%--------------------------------------------------------------------------
% Read the number of mesh nodes, nodal coordinates and boundary conditions.
%--------------------------------------------------------------------------
function [GEOM,BC,FEM] = innodes(GEOM,fid,FEM)
GEOM.npoin = fscanf(fid,'%d',1);     
format = ['%d %d ' repmat('%g ',1,GEOM.ndime)];
info = (fscanf(fid,format,[2+GEOM.ndime,GEOM.npoin]))';
nodes = info(:,1);
BC.icode = info(nodes,2);          
GEOM.x = info(nodes,3:end)'; 

for i = 1:FEM(1).n_elet_type
    FEM(i).mesh.n_dofs = GEOM.npoin*GEOM.ndime;           
    FEM(i).mesh.dof_nodes = reshape(1:GEOM.npoin*GEOM.ndime,GEOM.ndime,[]);  
end


end 
