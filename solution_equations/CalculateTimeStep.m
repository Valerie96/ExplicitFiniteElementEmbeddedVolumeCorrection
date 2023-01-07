%--------------------------------------------------------------------------
% Check for equilibrium convergence.
%--------------------------------------------------------------------------
function [dt] = CalculateTimeStep(FEM,GEOM,MAT,DAMPING)

dtMax = 1e20;
% Main element loop.
%--------------------------------------------------------------------------
for jmesh=1:length(FEM)
    switch FEM(jmesh).mesh.element_type 
        case'hexa8'
            
            for ielement=1:FEM(jmesh).mesh.nelem
                mu = MAT(jmesh).props(2);
                lambda = MAT(jmesh).props(3);
                rho = MAT(jmesh).props(1);
                %|-/
                b1 = DAMPING.b1;
                b2 = DAMPING.b2;
                eps_dot = GEOM.VolRate(ielement);
            
                %Longitudinal Bulk wave speed
                ce = sqrt((lambda + 2*mu)/rho);
                
%                 [le, max_le] = calc_element_size(FEM(jmesh),GEOM,ielement);
                global_nodes    = FEM(jmesh).mesh.connectivity(:,ielement);
                le = calc_min_element_size(GEOM.x(:,global_nodes),FEM(jmesh).mesh.element_type);
                dt_ielt = le/ce;
                
                %Add effect of damping 
                zeta = b1-(b2^2)*dt_ielt*min(0, eps_dot);
                dt_ielt = dt_ielt*(sqrt(1+zeta^2)-zeta);
                
                if(dt_ielt < dtMax)
                    dt = dt_ielt;
                    dtMax = dt;
                end
                
            end
           

        case 'truss2'

            for ielement=1:FEM(jmesh).mesh.nelem
                mu = MAT(jmesh).props(2);
                lambda = MAT(jmesh).props(3);
                rho = MAT(jmesh).props(1);
            
                %Longitudinal Bulk wave speed
                ce = sqrt((lambda + 2*mu)/rho);
                
%                 [le, max_le] = calc_element_size(FEM(jmesh),GEOM,ielement);
                global_nodes    = FEM(jmesh).mesh.connectivity(:,ielement);
                le = calc_min_element_size(GEOM.x(:,global_nodes),FEM(jmesh).mesh.element_type);
                dt_ielt = le/ce;
                
                if(dt_ielt < dtMax)
                    dt = dt_ielt;
                    dtMax = dt;
                end
                
            end
            
    end
end

end



