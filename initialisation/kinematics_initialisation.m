%--------------------------------------------------------------------------    
% Initialisation of kinematics. 
%--------------------------------------------------------------------------
function [KINEMATICS,INITIAL_KINEMATICS] = kinematics_initialisation(GEOM,FEM,QUADRATURE)

%Initialized the stucture arrays
KINEMATICS.Ib=0;
INITIAL_KINEMATICS.DN_X=zeros(3);

KINEMATICS = repmat(KINEMATICS,FEM(1).n_elet_type,1);
INITIAL_KINEMATICS = repmat(INITIAL_KINEMATICS,FEM(1).n_elet_type,1);

for i = 1:FEM(1).n_elet_type

    if strcmp(FEM(i).mesh.element_type,'truss2') || strcmp(FEM(i).mesh.element_type,'truss2_these_are_the_true_dim')
%          INITIAL_KINEMATICS(i) = [];
%          KINEMATICS(i) = [];
  
    else
        % Spatial gradient of the shape functions.
        KINEMATICS(i).DN_x   = zeros(GEOM.ndime,FEM(i).mesh.n_nodes_elem,QUADRATURE(i).element.ngauss);  
        % Jacobian of the mapping between spatial and isoparametric domains.
        KINEMATICS(i).Jx_chi = zeros(QUADRATURE(i).element.ngauss,1);  
        % Deformation gradient.
        KINEMATICS(i).F      = zeros(GEOM.ndime,GEOM.ndime,QUADRATURE(i).element.ngauss);                            
        % Left Cauchy-Green strain tensor (b).
        KINEMATICS(i).b      = zeros(GEOM.ndime,GEOM.ndime,QUADRATURE(i).element.ngauss);   
        % First invariant of b.
        KINEMATICS(i).Ib     = zeros(QUADRATURE(i).element.ngauss,1);      
        % Principal stretches.
        KINEMATICS(i).lambda = zeros(GEOM.ndime,QUADRATURE(i).element.ngauss);                
        % Spatial principal directions.
        KINEMATICS(i).n      = zeros(GEOM.ndime,GEOM.ndime,QUADRATURE(i).element.ngauss);
        % Previous Deformation Gradient Needs J to initially be 1
        KINEMATICS(i).J      = ones(QUADRATURE(i).element.ngauss,1);

        %-------------------------------------------------------------
        %Initilize Kinematics stucture to hold transition from natural
        %to initial coorinates space
        INITIAL_KINEMATICS(i).DN_X   = cell(FEM(i).mesh.nelem,1);  
        INITIAL_KINEMATICS(i).JX_chi = ones(FEM(i).mesh.nelem,QUADRATURE(i).element.ngauss,1);  

        INITIAL_KINEMATICS(i) = isoparametric_gradients(GEOM,FEM(i).mesh,FEM(i).interpolation.element.DN_chi,QUADRATURE(i).element,INITIAL_KINEMATICS(i));
    
    end
end
        
