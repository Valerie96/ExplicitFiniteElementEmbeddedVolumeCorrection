%--------------------------------------------------------------------------
% Computes J and the derivative of the shape functions with
% respect to the undeformed coordinates and its determinant (Jacobian of the 
% mapping between undeformed and isoparametric configurations) given the 
% initial coordinates of an element and the isoparametric derivative of
% the shape functions.
%--------------------------------------------------------------------------
function INITIAL_KINEMATICS = isoparametric_gradients(GEOM,mesh,DN_chi,QUADRATURE,INITIAL_KINEMATICS)
                                        
    for ielement=1:mesh.nelem
        INITIAL_KINEMATICS.DN_X{ielement,1}=zeros(GEOM.ndime,mesh.n_nodes_elem,QUADRATURE.ngauss);
        global_nodes    = mesh.connectivity(:,ielement);   
        Xlocal          = GEOM.x(:,global_nodes); 
        for igauss=1:QUADRATURE.ngauss
            %----------------------------------------------------------------------
            % Derivative of shape functions with respect to ...
            % - initial coordinates.
            %----------------------------------------------------------------------
            DX_chi = Xlocal*DN_chi(:,:,igauss)';
            DN_X   = DX_chi'\DN_chi(:,:,igauss);
            JX_chi = abs(det(DX_chi));
    
            %----------------------------------------------------------------------
            % Storage of variables.
            %----------------------------------------------------------------------
            INITIAL_KINEMATICS.DN_X{ielement,1}(:,:,igauss) = DN_X;  
            INITIAL_KINEMATICS.JX_chi(ielement,igauss)   = JX_chi;   
    
        end
    end
