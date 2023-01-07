%--------------------------------------------------------------------------
% Read and construct element type data. 
%--------------------------------------------------------------------------
function [FEM,GEOM,QUADRATURE] = elinfo(fid)                
FEM.n_elet_type                           = fscanf(fid,'%d',1);
QUADRATURE.element                        = [];
strtrim(fgetl(fid));

% create arrays of structures to hold data for multiple element types
FEM = repmat(FEM,FEM.n_elet_type,1);
QUADRATURE = repmat(QUADRATURE,FEM(1).n_elet_type,1);

for i = 1:FEM(1).n_elet_type
    FEM(i).mesh.element_type = strtrim(fgetl(fid));
    
    switch FEM(i).mesh.element_type
        case 'truss2'
             GEOM.ndime                               = 3;
             FEM(i).mesh.n_nodes_elem                 = 2;        
             QUADRATURE(i).element.polynomial_degree  = 1; 
             QUADRATURE(i).boundary                   = [];         
             FEM(i).mesh.n_face_nodes_elem            = 0;
             FEM(i).mesh.n_face_dofs_elem             = 0; 
        case 'tria3'
             GEOM.ndime                               = 2;
             FEM(i).mesh.n_nodes_elem                 = 3;
             FEM(i).mesh.n_face_nodes_elem            = 2;
             QUADRATURE(i).element.polynomial_degree  = 1; 
             QUADRATURE(i).boundary.polynomial_degree = 1;          
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
        case 'tria6'
             GEOM.ndime                               = 2;
             FEM(i).mesh.n_nodes_elem                 = 6;
             FEM(i).mesh.n_face_nodes_elem            = 3;
             QUADRATURE(i).element.polynomial_degree  = 2; 
             QUADRATURE(i).boundary.polynomial_degree = 2; 
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
        case 'quad4'
             GEOM.ndime                               = 2;
             FEM(i).mesh.n_nodes_elem                 = 4;
             FEM(i).mesh.n_face_nodes_elem            = 2;
             QUADRATURE(i).element.polynomial_degree  = 1; 
             QUADRATURE(i).boundary.polynomial_degree = 1; 
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
        case 'tetr4'
             GEOM.ndime                               = 3;
             FEM(i).mesh.n_nodes_elem                 = 4;
             FEM(i).mesh.n_face_nodes_elem            = 3;
             QUADRATURE(i).element.polynomial_degree  = 1; 
             QUADRATURE(i).boundary.polynomial_degree = 1; 
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
        case 'tetr10'
             GEOM.ndime                               = 3;
             FEM(i).mesh.n_nodes_elem                 = 10;
             FEM(i).mesh.n_face_nodes_elem            = 6;
             QUADRATURE(i).element.polynomial_degree  = 2; 
             QUADRATURE(i).boundary.polynomial_degree = 2; 
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
        case 'hexa8'
             GEOM.ndime                               = 3;
             FEM(i).mesh.n_nodes_elem                 = 8;
             FEM(i).mesh.n_face_nodes_elem            = 4;
             QUADRATURE(i).element.polynomial_degree  = 1; 
             QUADRATURE(i).boundary.polynomial_degree = 1; 
             FEM(i).mesh.n_face_dofs_elem             = ...
             FEM(i).mesh.n_face_nodes_elem*GEOM.ndime; 
    end    
    
    FEM(i).mesh.n_dofs_elem = FEM(i).mesh.n_nodes_elem*GEOM.ndime;   

end

end

