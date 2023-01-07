%----------------------------------------------------------------------
% Obtain fixed and free degree of freedom numbers (dofs).
%--------------------------------------------------------------------------
function BC = find_fixed_free_dofs(GEOM,FEM,BC)
%--------------------------------------------------------------------------
% Based on the boundary codes (BC.icode), determine the (free) degrees 
% of freedom (freedof) and the fixed degress of freedom (fixdof).
%--------------------------------------------------------------------------
dim = GEOM.ndime;
BC.freedof = (1:FEM.mesh.n_dofs)';
BC.fixdof = zeros(FEM.mesh.n_dofs,1);

BC.tienodes = zeros(GEOM.npoin,1);
BC.tiedof = zeros(FEM.mesh.n_dofs,1);
BC.hostdof = (1:FEM.mesh.n_dofs)';

% Sort fixed and free degrees of freedom
for inode=1:GEOM.npoin
    BC.fixdof((inode-1)*dim+(1:dim)) = boundary_codes(BC.icode(inode),dim);
end
BC.fixdof_truefalse = BC.fixdof;
BC.freedof_truefalse = ~BC.fixdof;

BC.fixdof = BC.fixdof.*BC.freedof;
BC.fixdof = BC.fixdof(BC.fixdof>0);

% Sort which dof are host nodes and which are embedded
for inode=1:GEOM.npoin
    if BC.icode(inode)==8
        BC.tiedof((inode-1)*dim+(1:dim)) = [1;1;1];
        BC.tienodes(inode)=inode;
        BC.hostdof((inode-1)*dim+(1:dim))=[0;0;0];
    end
end
BC.tiedof = BC.tiedof.*(1:FEM.mesh.n_dofs)';
BC.tiedof = BC.tiedof(BC.tiedof>0);
BC.tienodes = BC.tienodes(BC.tienodes>0);
BC.hostdof = BC.hostdof(BC.hostdof>0);
x = [BC.fixdof' BC.tiedof'];
BC.freedof(x')=[];


end



