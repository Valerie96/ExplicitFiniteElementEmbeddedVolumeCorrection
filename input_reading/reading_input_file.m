%--------------------------------------------------------------------------
% Read input data.
%--------------------------------------------------------------------------
function [FEM,GEOM,QUADRATURE,BC,MAT,LOAD,PRO,CON,GLOBAL,EmbedElt,...
    VolumeCorrect] = reading_input_file(PRO,CON,fid)
%--------------------------------------------------------------------------
% Problem title.   
%--------------------------------------------------------------------------
PRO.title = strtrim(fgets(fid));
%--------------------------------------------------------------------------
%Options: Embedded Elements, Volume Correction 
%--------------------------------------------------------------------------
text = fgetl(fid);
EmbedElt = sscanf(text, "Embedded Elements	%u"); 
text = fgetl(fid);
VolumeCorrect = sscanf(text, "VolumeCorrection	%u");    

%Read node and element for individual output
format = ['Output: Node %d Element %d'];
info = sscanf(fgetl(fid),format,2);
PRO.outnode  = info(1);
PRO.outelt   = info(2);
%--------------------------------------------------------------------------
% Element type.    
%--------------------------------------------------------------------------
[FEM,GEOM,QUADRATURE] = elinfo(fid);    
%--------------------------------------------------------------------------
% Obtain quadrature rules, isoparametric shape functions and their  
% derivatives for the internal and boundary elements.
%--------------------------------------------------------------------------
for i = 1:FEM(1).n_elet_type
      QUADRATURE(i).element = element_quadrature_rules(FEM(i).mesh.element_type);
      QUADRATURE(i).boundary = edge_quadrature_rules(FEM(i).mesh.element_type);

      FEM(i).interpolation = [];
      FEM(i) = shape_functions_iso_derivs(QUADRATURE(i),FEM(i),GEOM.ndime);
end
%--------------------------------------------------------------------------
% Read the number of mesh nodes, nodal coordinates and boundary conditions.  
%--------------------------------------------------------------------------
[GEOM,BC,FEM] = innodes(GEOM,fid,FEM);
%--------------------------------------------------------------------------
% Read the number of elements, element connectivity and material number.
%--------------------------------------------------------------------------
GEOM.total_n_elets = 0;
MAT.matno=zeros(2,1); MAT.props=zeros(2,1);MAT.nmats=0;MAT.n_nearly_incompressible=0;
MAT.matyp=0;
% create an array of MAT structures to hold data for multiple element types
MAT = repmat(MAT,FEM(1).n_elet_type,1);
for i = 1:FEM(1).n_elet_type
    [FEM(i),MAT(i)] = inelems(FEM(i),MAT(i) ,fid);
    GEOM.total_n_elets = GEOM.total_n_elets + FEM(i).mesh.nelem;
end
%--------------------------------------------------------------------------
% Obtain fixed and free degree of freedom numbers (dofs).
%--------------------------------------------------------------------------
BC = find_fixed_free_dofs(GEOM,FEM(1),BC);
%--------------------------------------------------------------------------
% Read the number of materials and material properties.  
%--------------------------------------------------------------------------
for i = 1:FEM(1).n_elet_type
    MAT(i) = matprop(MAT(i),FEM(i),fid);                          
end
%--------------------------------------------------------------------------
% Read nodal point loads, prescribed displacements, surface pressure loads
% and gravity (details in textbook).
%--------------------------------------------------------------------------
[LOAD,BC,FEM(1),GLOBAL,CON] = inloads(GEOM,FEM(1),BC,CON,fid);
%--------------------------------------------------------------------------

fclose('all'); 
