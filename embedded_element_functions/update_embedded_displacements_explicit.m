%--------------------------------------------------------------------------
%  Update coodinates of embedded element displacements.
%--------------------------------------------------------------------------       
function [x,v,a]       = update_embedded_displacements_explicit(tiedof,...
                    tienodes,FEM,GEOM,v,a)
mesh = FEM(2).mesh;
%Enforce Embedded Element constraint

    %-------------------------------------------------------------------------- 
    %Use reference coordinates of host elt to establish embedded element
    %natural coordinates. The updated location is found using current nodal
    %coordinatess
    %-------------------------------------------------------------------------- 
    dim = GEOM.ndime;
    dof = mesh.dof_nodes;
    x = GEOM.x;
    Ze = GEOM.embedded.Embed_Zeta;
    TieXUpdate = zeros(mesh.n_dofs,1);
    TieVUpdate = zeros(mesh.n_dofs,1);
    TieAUpdate = zeros(mesh.n_dofs,1);

    %-------------------------------------------------------------------------- 
    %Loop through embedded nodes, m
    %Disp at a tie node is = displacment of that point in the host elet
    %-------------------------------------------------------------------------- 
    for i=1:length(tienodes)
        m=tienodes(i);

        %-------------------------------------------------------------------------- 
        %Get the current coordinates of the host element
        host = GEOM.embedded.NodeHost(m);        %host element number
        host_nn=FEM(1).mesh.connectivity(:,host); %nodes of host element
        host_xn = x(:,host_nn);            %nodal coordinates of embedded elet in host
        host_dof = dof(:, host_nn);
        %-------------------------------------------------------------------------- 
        %Get shape function values at embedded node location
        N=shape_function_values_at(Ze(:,m),'hexa8')';
        %-------------------------------------------------------------------------- 
        %Calculate interpolated displacment values at embedded node location
         XeUpdate=N*host_xn';
        %-------------------------------------------------------------------------- 
        %Calculate interpolated velocity values
         V_Update=N*v(host_dof)';
        %-------------------------------------------------------------------------- 
        %Calculate interpolated acceleration values
         A_Update=N*a(host_dof)';
        %-------------------------------------------------------------------------- 
        %Round machine zero results to actual zero
        for j=1:length(XeUpdate)
            if abs(XeUpdate(j)) < 1E-12
                XeUpdate(j)=0;
            end
            if abs(V_Update(j)) < 1E-12
                V_Update(j)=0;
            end
            if abs(A_Update(j)) < 1E-12
                A_Update(j)=0;
            end
        end

        %Fill TieUpdate with new xyz locations
        TieXUpdate((m-1)*dim+(1:dim)) = XeUpdate;
        TieVUpdate((m-1)*dim+(1:dim)) = V_Update;
        TieAUpdate((m-1)*dim+(1:dim)) = A_Update;
    end

    x(tiedof) = TieXUpdate(tiedof);
    v(tiedof) = TieVUpdate(tiedof);
    a(tiedof) = TieAUpdate(tiedof);

end