function GEOM = inverse_mapping(GEOM,FEM,tienodes)
% Calculates the host element coordinates of embedded elements and creates
% NodeHost: list of nodes and their associated host, and ElementHost: list
% of elements and their associated host. 


%-----------------------------------------------------------------------
%!!!!!!!!!!!
%Assume all elements of the first type are hosts and all of the
%second type are embedded
%!!!!!!!!!!!
%-----------------------------------------------------------------------
    e_nodes = tienodes(:);
    h_elts=FEM(1).mesh.host; 
    e_elts=FEM(2).mesh.embedded;


    NodeHost = zeros(GEOM.npoin,1);
    ElementHost = zeros(FEM(2).mesh.nelem,4);
    HostTotals = zeros(FEM(1).mesh.nelem,2);
    HostsElements = cell(FEM(1).mesh.nelem,1);
    Zeta = zeros(3,GEOM.npoin);
    
    %Loop throgh host elements
    for j = 1:length(h_elts)
        h = h_elts(j);
        h_connectivity = FEM(1).mesh.connectivity(:,h);
        x_h = GEOM.x0(:,h_connectivity);  %Host node global coordinates
        HostsElements{j}=zeros(FEM(2).mesh.nelem,1);

        %Loop over all embedded nodes
        for i = 1:length(e_nodes)
            ne = e_nodes(i);

            %Check if we already did this node
            if NodeHost(ne) == 0
                x_ne = GEOM.x0(:,ne);
                inel = point_in_hexahedron(x_ne', x_h);

                %Check if it's in this host elt
                if inel
                    NodeHost(ne) = h;
                    HostTotals(h,1) = HostTotals(h,1) + 1;

                    %Find natrual coordinates of the node
                    Zeta(:,ne) = find_natural_coords(x_ne, x_h, FEM(1).mesh.element_type);
                end
            end
        end
    end

    %Assign hosts to embedded elements
    for i=1:length(e_elts)
        e = e_elts(i); 
        e_connectivity = FEM(2).mesh.connectivity(:,e);
        n1 = e_connectivity(1); %Choose first node
        host1 = NodeHost(n1);
        ElementHost(e,1) = host1;
        n2 = e_connectivity(2); %Choose second node
        host2 = NodeHost(n2);
        ElementHost(e,2) = host2;
        
        if host1==host2 && host1~=0
            HostTotals(host1,2) = HostTotals(host1,2) + 1;
            ElementHost(e,3) = 0.5;
            ElementHost(e,4) = 0.5;
            HostsElements{host1}(e)=e;
        else
            p = find_intersection(GEOM.x0,FEM,e_connectivity,host1,host2);
            
            %Persentage of truss in each host
            %truss vector
                truss = GEOM.x0(:,e_connectivity(2)) - GEOM.x0(:,e_connectivity(1));
            %vector from p to truss node 1
                truss_h1 = p - GEOM.x0(:,e_connectivity(1));
            %percent of total host length between point p and node 1
                percent1 = norm(truss_h1)/norm(truss);
            ElementHost(e,3) = percent1;
            ElementHost(e,4) = 1-percent1;
            
            HostTotals(host1,2) = HostTotals(host1,2) + 1;
            HostTotals(host2,2) = HostTotals(host2,2) + 1;
            HostsElements{host1}(e)=e;
            HostsElements{host2}(e)=e;
            
        end

    end

    GEOM.embedded.NodeHost = NodeHost;
    GEOM.embedded.ElementHost = ElementHost;
    GEOM.embedded.HostTotals = HostTotals;
    GEOM.embedded.Embed_Zeta = Zeta;

    for i=1:FEM(1).mesh.nelem
        HostsElements{i} = HostsElements{i}(HostsElements{i}>0);
    end
    GEOM.embedded.HostsElements=HostsElements;
end