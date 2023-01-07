%--------------------------------------------------------------------------
% Computes the element vector of global internal forces 
%--------------------------------------------------------------------------
function [T_internal,geomJn_1,VolRate,f_damp,temp_STRESS] = ...
          InternalForce_explicit(ielement,FEM,xlocal,x0local,...
          Ve,QUAD,properties,CONSTANT,GEOM,...
          PLAST,matyp,KINEMATICS,DN_X,MAT,DAMPING,STRESS,VolumeCorrect,dt)
      
dim=GEOM.ndime;

% step 2.II       
T_internal = zeros(FEM(1).mesh.n_dofs_elem,1);
f_damp     = zeros(FEM(1).mesh.n_dofs_elem,1);
temp_STRESS = struct('Cauchy',zeros(1,6,8),'LE',zeros(1,6,8));

%--------------------------------------------------------------------------
% Computes initial and current gradients of shape functions and various 
% strain measures at all the Gauss points of the element.
%--------------------------------------------------------------------------
KINEMATICS(1) = gradients(xlocal,x0local,FEM(1).interpolation.element.DN_chi,...
             QUAD,KINEMATICS(1),DN_X);


Jn_1=GEOM.Jn_1(ielement);
J=KINEMATICS(1).J(1);
eps_dot = (J-Jn_1)/dt;
b1=DAMPING.b1;
b2=DAMPING.b2;

%--------------------------------------------------------------------------
% Computes element mean dilatation kinematics, pressure and bulk modulus. 
%--------------------------------------------------------------------------
switch matyp
     case {5,7,17}
          [pressure,kappa_bar,DN_x_mean,ve] = ...
           mean_dilatation_pressure(FEM(1).mesh.n_nodes_elem,dim,matyp,properties,Ve,...
                                    QUAD,KINEMATICS(1));
     otherwise
          pressure = 0;
end
%--------------------------------------------------------------------------
% Gauss quadrature integration loop.
%--------------------------------------------------------------------------
for igauss=1:QUAD.ngauss
    %----------------------------------------------------------------------
    % Extract kinematics at the particular Gauss point.
    %----------------------------------------------------------------------
    kinematics_gauss = kinematics_gauss_point(KINEMATICS(1),igauss);     
    %----------------------------------------------------------------------
    % Obtain stresses (for incompressible or nearly incompressible, 
    % only deviatoric component) and internal variables in plasticity.
    %----------------------------------------------------------------------    
    [Cauchy,~,~] = Cauchy_type_selection(kinematics_gauss,properties,...
                                          CONSTANT,dim,matyp,PLAST,igauss);
    %----------------------------------------------------------------------
    % Elasticity tensor is not used in explicit analysis
    %----------------------------------------------------------------------   
    c = 0;
    %----------------------------------------------------------------------
    % Add pressure contribution to stresses and elasticity tensor.
    %----------------------------------------------------------------------    
    [Cauchy,~] = mean_dilatation_pressure_addition(Cauchy,c,CONSTANT,pressure,matyp); 
    
    %----------------------------------------------------------------------
    % Calculate bulk viscosity damping
    global_nodes=FEM(1).mesh.connectivity(:,ielement);
    le=calc_min_element_size(GEOM.x(:,global_nodes),FEM(1).mesh.element_type);
    rho=properties(1); mu=properties(2); lambda=properties(3);
    Cd=sqrt((lambda + 2*mu)/rho);
    
    p1 = rho*b1*le*Cd*eps_dot*CONSTANT.I;
    p2 = rho*(b2*le)^2*abs(eps_dot)*min(0,eps_dot)*CONSTANT.I;
        
    %----------------------------------------------------------------------
    % Compute numerical integration multipliers.
    %----------------------------------------------------------------------
    JW = kinematics_gauss.Jx_chi*QUAD.W(igauss)*...
         thickness_plane_stress(properties,kinematics_gauss.J,matyp);
    %----------------------------------------------------------------------
    % Compute equivalent (internal) force vector.
    %----------------------------------------------------------------------
    T = (Cauchy+p1+p2)*kinematics_gauss.DN_x;
    T_internal = T_internal + T(:)*JW;

    fd = (p1+p2)*kinematics_gauss.DN_x; 
    f_damp = f_damp + fd(:)*JW;

    %Save element stress and strain
       components = [1;2;3;5;6;9];
       %-----------------------------------------------------------------------
       %I want log strain too
       %-----------------------------------------------------------------------
       lam = kinematics_gauss.lambda;
       LE = zeros(dim,dim);
       for j=1:dim
           LE = LE + log(lam(j))*kinematics_gauss.n(:,j)*kinematics_gauss.n(:,j)'; 
           LE(~isfinite(LE)) = 0;
       end
       temp_STRESS.Cauchy(1,:,igauss) = Cauchy(components);
       temp_STRESS.LE(1,:,igauss) = LE(components);

        
end

    %Update previous Jacobian and element strain rate
    %Assuming that J and eps_dot are the same for all the element Gauss Pts
    %and I have now confirmed this
    geomJn_1=J;
    VolRate = eps_dot;

    
%--------------------------------------------------------------------------
% Compute conttribution (and extract relevant information for subsequent
% assembly) of the mean dilatation term (Kk) of the stiffness matrix.
%--------------------------------------------------------------------------
% switch matyp
%     case {5,7,17}         
%          [indexi,indexj,global_stiffness,...
%           counter] = mean_dilatation_volumetric_matrix(FEM,dim,...
%           element_connectivity,DN_x_mean,counter,indexi,indexj,...
%           global_stiffness,kappa_bar,ve);
% end 


% Embedded Elt, Internal force modification, if this element has any
% embedded elements
if GEOM.embedded.HostTotals(ielement,2) > 0
    for eelt=GEOM.embedded.HostsElements{ielement}'
        node_flag = [0 0];
        if GEOM.embedded.ElementHost(eelt,1) == ielement
            node_flag(1) = 1;
        end
        if GEOM.embedded.ElementHost(eelt,2) == ielement
            node_flag(2) = 1;
        end
        if node_flag(1) + node_flag(2) >= 1;
            material_number   = MAT(2).matno(eelt);
            properties_e      = MAT(2).props(:,material_number);
            T_internal = TrussCorrectedInternalForce_explicit_from_mem(ielement,...
                           T_internal,FEM(2).mesh.connectivity,FEM(1).mesh.element_type,GEOM,PLAST,STRESS(2),...
                           properties,properties_e,DAMPING,VolumeCorrect,eelt,node_flag);
        end
    end
end


%--------------------------------------------------------------------------
% Store internal variables.
%--------------------------------------------------------------------------
% PLAST_element = PLAST;
PLAST_element=[];
end

