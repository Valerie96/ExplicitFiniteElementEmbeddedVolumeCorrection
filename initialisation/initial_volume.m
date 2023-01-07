%--------------------------------------------------------------------------    
% Calculate initial volume for data checking. 
%--------------------------------------------------------------------------
function GEOM = initial_volume(FEM,GEOM,QUADRATURE,MAT,INITIAL_KINEMATICS)

%Determine the largest number of elements 
max_elet_type = 0;
for i=1:FEM(1).n_elet_type
    if FEM(i).mesh.nelem > max_elet_type
        max_elet_type = FEM(i).mesh.nelem;
    end
end

Ve      = zeros(max_elet_type,FEM(1).n_elet_type); 
V_total = 0;
M_total = 0;


for i = 1:FEM(1).n_elet_type
    for ielement=1:FEM(i).mesh.nelem
        global_nodes    = FEM(i).mesh.connectivity(:,ielement);   
        xlocal          = GEOM.x(:,global_nodes);              
        material_number = MAT(i).matno(ielement);                 
        matyp           = MAT(i).matyp(material_number);          
        properties      = MAT(i).props(:,material_number); 
        
        switch FEM(i).mesh.element_type
            case 'truss2'        
                 area         = properties(4);  
                 L            = norm(xlocal(:,2) - xlocal(:,1));    
                 Ve(ielement,i) = area*L;                             
            otherwise
                 %-------------------------------------------------------------
                 % Compute gradients with respect to isoparametric coordinates.
                 %-------------------------------------------------------------
                 for igauss = 1:QUADRATURE(i).element.ngauss
                     %---------------------------------------------------------
                     % Computes the thickness in the deformed configuration for
                     % plane stress problems.
                     %---------------------------------------------------------
                     thickness_factor  = thickness_plane_stress(properties,1,matyp);%Initial Jacobian between current and deformed is 1 because there is no deformation yet            
                     JW = INITIAL_KINEMATICS(i).JX_chi(igauss)*QUADRATURE(i).element.W(igauss)*thickness_factor; 
                     %---------------------------------------------------------
                     % Compute volume of an element of the mesh.
                     %---------------------------------------------------------
                     Ve(ielement,i) = Ve(ielement,i) + JW;
                 end
        end
        %----------------------------------------------------------------------
        % Total volume of the mesh. 
        V_total = V_total + Ve(ielement,i);
        %Total Mass of mesh (for fun)
        M_total = M_total + Ve(ielement,i)*properties(1);
    end
end
%--------------------------------------------------------------------------
% Save information in data structure.
%--------------------------------------------------------------------------
GEOM.Ve      = Ve;
GEOM.V_total = V_total;
fprintf('Total mesh volume is: %15.5f \n', GEOM.V_total)
fprintf('Total mesh mass (p*V) is: %15.5f \n', M_total)



