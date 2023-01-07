%Calculates the embedded element effect on host element during
%InternalForce_explicit


function [T_internal] = TrussCorrectedInternalForce_explicit_from_mem(ielement,...
          T_internal,connectivity2,host_element_type,GEOM,PLAST,STRESS,...
          properties_h,properties_e,DAMPING,VolumeCorrect,eelt,node_flag)

dim=GEOM.ndime;

% *** For now we are assuming that host elements are of FEM type 1 and
% embedded are FEM type 2
    %----------------------------------------------------------------------
    % GATHER material properties of the host and embedded elt
    %----------------------------------------------------------------------   
%     material_number   = MAT(1).matno(ielement);     
%     matyp_h           = MAT(1).matyp(material_number);        
%     properties_h      = MAT(1).props(:,material_number);                       
%     Ve_h              = GEOM.Ve(ielement);  
    
%     material_number   = MAT(2).matno(eelt);     
%     matyp_e           = MAT(2).matyp(material_number);        
%     properties_e      = MAT(2).props(:,material_number);
%     Ve_e              = GEOM.Ve(eelt);  
    
    %Get embedded element information
    e_connectivity = connectivity2(:,eelt);%FEM(2).mesh.connectivity(:,eelt);
    x_e = GEOM.x0(:,e_connectivity);
    xelocal  = GEOM.x(:,e_connectivity);                     
    e_nodes_zeta = GEOM.embedded.Embed_Zeta(:,e_connectivity);
    
    %Node flag indicates which nodes are actually inside the host for the
    %case where one truss element could span multiple hosts (1 for in, 0
    %for out)
    percent1 = node_flag(1);
    percent2 = node_flag(2);
     
    %--------------------------------------------------------------------------
    % Get Correction force (Embedded Truss internal forces using host element
    % material properties)
    %--------------------------------------------------------------------------    
    %Assuming both elements are neo hookean materials. But the truss model
    %requires E and nu instead of lambda and mu
    lam_h=properties_h(3);
    mu_h=properties_h(2);
    K_h=lam_h+(2*mu_h/3);
    nu_h   = (3*K_h-2*mu_h)/(2*(3*K_h+mu_h));
    E_h    = 9*K_h*mu_h/(3*K_h+mu_h);
    
    properties_eh = properties_e; %eh is the same as e, except that nu and E are replaced bu nu and E of the host element
    properties_eh(1) = properties_h(1);
    properties_eh(2) = E_h;
    properties_eh(3) = nu_h;

    TE=STRESS.InternalForce(eelt,:);
    [TC,~,~,~,~,~,~] = element_force_truss(properties_eh,xelocal,x_e,PLAST,GEOM,DAMPING,1);

    %----------------------------------------------------------------------
    % Get embeddded element internal force and convert to force on host
    % nodes
    %---------------------------------------------------------------------- 
    % Shape functions at embedded element node in their respective hosts
    % (not nessisarily ielement)
    N_node1 = shape_function_values_at(e_nodes_zeta(:,1), host_element_type);
    N_node2 = shape_function_values_at(e_nodes_zeta(:,2), host_element_type);
    
    %Force from embedded nodes distribted over host nodes
    T_e1 = zeros(dim*8, 1);
    T_e2 = zeros(dim*8, 1);
    %Do the same for the correction force
    T_C1 = zeros(dim*8, 1);
    T_C2 = zeros(dim*8, 1);
    for i = 1:3:24
       T_e1(i:i+2) = TE(1:3)*N_node1((i-1)/3 + 1)*percent1*node_flag(1); 
       T_e2(i:i+2) = TE(4:6)*N_node2((i-1)/3 + 1)*percent2*node_flag(2); 

       T_C1(i:i+2) = TC(1:3)*N_node1((i-1)/3 + 1)*percent1*node_flag(1); 
       T_C2(i:i+2) = TC(4:6)*N_node2((i-1)/3 + 1)*percent2*node_flag(2);
    end
    T_e = T_e1 + T_e2; 
    T_C = T_C1 + T_C2;
      
    
    %----------------------------------------------------------------------
    % Compute equivilant (internal) force vector of the host element.
    %----------------------------------------------------------------------
    if VolumeCorrect
        T_internal = T_internal + (T_e - T_C);
    else
        T_internal = T_internal + (T_e);
    end
 
        
end