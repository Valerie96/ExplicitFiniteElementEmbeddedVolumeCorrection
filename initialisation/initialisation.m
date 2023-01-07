%----------------------------------------------------------------------
% Initialises kinematic variables and computes initial tangent matrix 
% and equivalent force vector, excluding pressure components.
%----------------------------------------------------------------------
function [GEOM,LOAD,GLOBAL,PLAST,KINEMATICS,INITIAL_KINEMATICS,STRESS] = ...
         initialisation(FEM,GEOM,QUADRATURE,MAT,LOAD,CON,GLOBAL,BC,EmbedElt,VolumeCorrect)
%--------------------------------------------------------------------------    
% Initialisation of internal variables for plasticity.
%--------------------------------------------------------------------------    
GEOM.element_num = zeros(3,GEOM.total_n_elets);
PLAST.a=[0];
PLAST = repmat(PLAST,FEM(1).n_elet_type,1);

% Initilize structure to store stress, strain, and truss internal force
STRESS.InternalForce = [0];
STRESS = repmat(STRESS,FEM(1).n_elet_type,1);

for i = 1:FEM(1).n_elet_type
check = (isempty((MAT(i).matyp(MAT(i).matyp==17)))*isempty((MAT(i).matyp(MAT(i).matyp==2))));

    if check  
       PLAST(i).a = [0];
       PLAST(i).b = [0];

       switch FEM(i).mesh.element_type
           case 'truss2'
                STRESS(i).Cauchy=zeros(FEM(i).mesh.nelem,1);
                STRESS(i).LE=zeros(FEM(i).mesh.nelem,1);
                STRESS(i).InternalForce=zeros(FEM(i).mesh.nelem,GEOM.ndime*2);
           otherwise
                STRESS(i).Cauchy=zeros(FEM(i).mesh.nelem,6,QUADRATURE(i).element.ngauss);
                STRESS(i).LE=zeros(FEM(i).mesh.nelem,6,QUADRATURE(i).element.ngauss);
       end

    else 
   
       switch FEM(i).mesh.element_type
           case 'truss2'
                PLAST(i).ep    = zeros(FEM(i).mesh.nelem,1);  
                PLAST(i).epbar = zeros(FEM(i).mesh.nelem,1); 
                STRESS(i).Cauchy=zeros(FEM(i).mesh.nelem,1);
                STRESS(i).LE=zeros(FEM(i).mesh.nelem,1);
                STRESS(i).InternalForce=zeros(FEM(i).mesh.nelem,GEOM.ndime*2);
           otherwise
                PLAST(i).epbar = zeros(QUADRATURE(i).element.ngauss,FEM(i).mesh.nelem,1);       
                PLAST(i).invCp = reshape(repmat(eye(GEOM.ndime),1,...
                              QUADRATURE(i).element.ngauss*FEM(i).mesh.nelem),...
                              GEOM.ndime,GEOM.ndime,QUADRATURE(i).element.ngauss,...
                              FEM(i).mesh.nelem); 
                STRESS(i).Cauchy=zeros(FEM(i).mesh.nelem,QUADRATURE(i).element.ngauss);
                STRESS(i).LE=zeros(FEM(i).mesh.nelem,QUADRATURE(i).element.ngauss);
       end
    end

%--------------------------------------------------------------------------
end
%Initialisation of kinematics. 
%--------------------------------------------------------------------------
[KINEMATICS,INITIAL_KINEMATICS] = kinematics_initialisation(GEOM,FEM,QUADRATURE);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------    
% Initialise undeformed geometry and initial residual and external forces. 
%--------------------------------------------------------------------------    
mesh_dof             = FEM(1).mesh.n_dofs;
GEOM.x0              = GEOM.x;
GLOBAL.Residual      = zeros(mesh_dof,1);
GLOBAL.external_load = zeros(mesh_dof,1);
GLOBAL.Reactions     = zeros(mesh_dof,1);

%--------------------------------------------------------------------------  
% Define velocity and accelerations
%--------------------------------------------------------------------------  
GLOBAL.velocities = zeros(mesh_dof,1);
GLOBAL.accelerations = zeros(mesh_dof,1);
GEOM.Jn_1 = ones(mesh_dof,1); %Jacobian of the previous step
GEOM.VolRate = zeros(mesh_dof,1);

%--------------------------------------------------------------------------       
% Calculate initial volume for data checking. 
% Additionally, essential for mean dilation algorithm.
%--------------------------------------------------------------------------
GEOM = initial_volume(FEM,GEOM,QUADRATURE,MAT,INITIAL_KINEMATICS);
%--------------------------------------------------------------------------    
% Compute the external force vector contribution due to gravity 
% (nominal value prior to load increment). 
%--------------------------------------------------------------------------    
if norm(LOAD.gravt)>0
   fprintf("gravity doesn't work yet\n");
   GLOBAL = gravity_vector_assembly(GEOM,FEM,QUADRATURE.element,LOAD,...
                                    MAT,GLOBAL,KINEMATICS);     
end
%--------------------------------------------------------------------------    
% Initialise external force vector contribution due to pressure
% (nominal value prior to load increment).
%--------------------------------------------------------------------------    
GLOBAL.nominal_pressure = zeros(FEM(1).mesh.n_dofs,1);
%--------------------------------------------------------------------------    
% Computes and assembles the initial external loads
%-------------------------------------------------------------------------- 
    GLOBAL.external_load = CON.xlamb*GLOBAL.nominal_external_load;
    GLOBAL.T_int         = zeros(FEM(1).mesh.n_dofs,1);
    GLOBAL.Residual      = GLOBAL.T_int - GLOBAL.external_load;

%--------------------------------------------------------------------------    
% Assemble the mass matrix and map embedded nodes/elements to hosts (if
% there are any)
%--------------------------------------------------------------------------     
    if (EmbedElt == 1)
        % map embedded nodes/elements to hosts
        GEOM             = inverse_mapping(GEOM,FEM,BC.tienodes);
        % assemble effective mass matrix for the host elements only 
        [GLOBAL]         = effective_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE,VolumeCorrect);
        GLOBAL.M         = GLOBAL.M(BC.hostdof(:,1),BC.hostdof(:,1));

    else
        % no mapping embedded nodes/elements to hosts
        GEOM.embedded.NodeHost    = zeros(GEOM.npoin,1);
        GEOM.embedded.ElementHost = zeros(FEM(1).mesh.nelem,8);
        GEOM.embedded.HostTotals  = zeros(FEM(1).mesh.nelem,2);
        GEOM.embedded.Embed_Zeta  = zeros(3, GEOM.npoin);
        % assemble mass matrix
        [GLOBAL] = lumped_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE);     
    end
       
end

