%--------------------------------------------------------------------------
% Obtain stresses (for incompressible or nearly incompressible, 
% only deviatoric component) and internal variables in plasticity.
%--------------------------------------------------------------------------
function [Cauchy,PLAST,PLAST_gauss] = Cauchy_type_selection(kinematics,...
          properties,cons,dim,matyp,PLAST,igauss)
PLAST_gauss = [];
switch matyp
    case {1,2} %This is just a sloppy change. Mat type 2 is for truss 
               %elements that really shouldn't come through this function
         Cauchy = stress1(kinematics,properties,cons);
    case 3
         Cauchy = stress3(kinematics,properties,cons);
end



