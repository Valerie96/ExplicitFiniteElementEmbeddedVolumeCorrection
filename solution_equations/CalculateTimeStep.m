%--------------------------------------------------------------------------
% Check for equilibrium convergence.
%--------------------------------------------------------------------------
function [dt] = CalculateTimeStep(FEM,GEOM,MAT,DAMPING)

dtMax = 1e20;

for jmesh=1:length(FEM)  
    matno_list      = MAT(jmesh).matno;     
    matyp_list      = MAT(jmesh).matyp;        
    props_list      = MAT(jmesh).props; 

            
            for ielement=1:FEM(jmesh).mesh.nelem 
                material_number = matno_list(ielement);     
                matyp           = matyp_list(material_number);        
                properties      = props_list(:,material_number);

                switch matyp
                    case  1
                        mu      = properties(2);
                        lambda  = properties(3);
                        rho     = properties(1);
                        %Longitudinal Bulk wave speed
                            ce = sqrt((lambda + 2*mu)/rho);
                    case 2
                        E       = properties(2);
                        %nu      = properties(3);
                        rho     = properties(1);
                        %Longitudinal Bulk wave speed
                           ce = sqrt(E/rho);
                    case 3
                        C10     = properties(2);
                        C01     = properties(3);
                        D1      = properties(4);
                        rho     = properties(1);

                        mu     = 2*(C10+C01);
                        K      = 2/D1;
                        lambda = K-2*mu/3;
                        %Longitudinal Bulk wave speed
                            ce = sqrt((lambda + 2*mu)/rho);
                    case 4 
                        rho   = properties(1);
                        alpha1     = properties(2);
                        mu1        = properties(3);
                        D1         = properties(4);
                        K=2/D1;
                        mu = mu1*alpha1/2;
                        lambda = K-2*mu/3;
                        %Longitudinal Bulk wave speed
                            ce = sqrt((lambda + 2*mu)/rho);
                    otherwise
                        fprintf("Material type not programed");
                end


                global_nodes    = FEM(jmesh).mesh.connectivity(:,ielement);
                le = calc_min_element_size(GEOM.x(:,global_nodes),FEM(jmesh).mesh.element_type);
                dt_ielt = le/ce;
                
                if strcmp(FEM(jmesh).mesh.element_type,'hexa8')
                    %Add effect of damping 
                    b1 = DAMPING.b1;
                    b2 = DAMPING.b2;
                    eps_dot = GEOM.VolRate(ielement);
    
                    zeta = b1-(b2^2)*dt_ielt*min(0, eps_dot);
                    dt_ielt = dt_ielt*(sqrt(1+zeta^2)-zeta);
                end
                
                if(dt_ielt < dtMax)
                    dt = dt_ielt;
                    dtMax = dt;
                end
                
            end
           
            
end


end



