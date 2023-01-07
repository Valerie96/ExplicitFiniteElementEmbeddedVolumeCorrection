function N = shape_function_values_at(Z, eltype)
%Returns shape function values when evaluated ad Z (vector of coordinates)

switch eltype
    case'hexa8'
    chi = Z(1); eta = Z(2); iota = Z(3);


     N1 = -((chi - 1)*(eta - 1)*(iota - 1))/ 8;
     N2 = ((chi + 1)*(eta - 1)*(iota - 1))/ 8;
     N3 = -((chi + 1)*(eta + 1)*(iota - 1))/ 8;
     N4 = ((chi - 1)*(eta + 1)*(iota - 1))/ 8;
     N5 = ((chi - 1)*(eta - 1)*(iota + 1))/ 8;
     N6 = -((chi + 1)*(eta - 1)*(iota + 1))/ 8;
     N7 = ((chi + 1)*(eta + 1)*(iota + 1))/8;
     N8 = -((chi - 1)*(eta + 1)*(iota + 1))/8;

     N = [N1; N2; N3; N4; N5; N6; N7; N8];
     
    case 'truss2'
        nu = Z;
        
        N1 = (1-nu)/2;
        N2 = (1+nu)/2;
        
        N = [N1; N2];
end