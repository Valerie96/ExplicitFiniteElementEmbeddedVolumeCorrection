%--------------------------------------------------------------------------
% Computes and assemble lumped mass matrix. Mass Correct gives the option
% to remove redundant mass from the system when using embedded elements.
%--------------------------------------------------------------------------
function [GLOBAL] = lumped_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE)   

%number of dimensions
ndims = GEOM.ndime; 
massSize = ndims*GEOM.npoin;
LumpedMass = zeros(massSize,massSize);


M_total = 0;
%Loop over all (host) elements
for ii=1:FEM(1).mesh.nelem
    ielement=ii;
    %----------------------------------------------------------------------
    % Temporary variables associated with a particular element.
    %----------------------------------------------------------------------
    n_nodes_elem    = FEM(1).mesh.n_nodes_elem;
    global_nodes    = FEM(1).mesh.connectivity(:,ielement);   
    material_number = MAT(1).matno(ielement);     
    matyp           = MAT(1).matyp(material_number);        
    properties      = MAT(1).props(:,material_number); % needs density included.
    xlocal          = GEOM.x(:,global_nodes);                     
    x0local         = GEOM.x0(:,global_nodes); 
    Ve              = GEOM.Ve(ielement,1);
   

    
    % assign density
    rho = properties(1);
    
    
    % quadrature locations
    QUADRATURE(1).element.Chi;
    % quadrature weights
    QUADRATURE(1).element.W;
   

    %b: calculate element mass
    me = rho * Ve;
        
        %d: calculate effective mass
        M_total = M_total + me;
        
        %e: divide mass equally among nodes
        Meff = (me/n_nodes_elem);
        
        %f: form diagonal mass matrix
        Me = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Meff;

        
        %g: transform and scatter scatter mass matrix to global 
        global_dof = FEM(1).mesh.dof_nodes(:,global_nodes);
        LumpedMass(global_dof,global_dof) = LumpedMass(global_dof,global_dof) + Me;
        
        
end % loop on elements

GLOBAL.M=LumpedMass;
GLOBAL.M_KE=LumpedMass;

M_total = 0 ;
for i = 1:massSize
    M_total = M_total + LumpedMass(i,i);
end

end

 


 
