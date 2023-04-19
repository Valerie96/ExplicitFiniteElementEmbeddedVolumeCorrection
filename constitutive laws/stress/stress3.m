%--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor for material type 3 (3D Mooney Rivlin).
%--------------------------------------------------------------------------
function Cauchy = stress3(kinematics,properties,cons)


C10             = properties(2);
C01             = properties(3);
D1              = properties(4);
J               = kinematics.J;
b               = kinematics.b;


%I bar formulation
K               = 2/D1;
I1              = kinematics.Ib;
I2              = (1/2)*(trace(b)^2-trace(b*b));
Cauchy          = J^(-5/3)*2*C10*(b-(1/3)*I1*cons.I) + J^(-7/3)*2*C01*(I1*b-b*b-(2/3)*I2*cons.I) + K*(J-1)*cons.I;

    
end