%--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor for material type 1.
%--------------------------------------------------------------------------
function Cauchy = stress1(kinematics,properties,cons)


mu              = properties(2);
lambda          = properties(3);
J               = kinematics.J;
b               = kinematics.b;
Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;

%I bar formulation
K               = (3*lambda+2*mu)/3;
I1              = trace(b);
Cauchy          = J^(-5/3)*mu*(b-(1/3)*I1*cons.I) + K*(J-1)*cons.I;

% % Linear Elastic
% F               = kinematics.F;
% % E          = (1/2)*(F'*F-eye(3));
% E  = (sqrtm(b)-eye(3));
% LogStrain = logm(sqrtm(b));
% E=LogStrain;
% E6=[E(1,1) E(2,2) E(3,3) 2*E(2,3) 2*E(1,3) 2*E(1,2)]';
% 
% lam=lambda;
% G=mu;
% C=[lam+2*G lam lam 0 0 0; lam lam+2*G lam 0 0 0; lam lam lam+2*G 0 0 0; 0 0 0 G 0 0; 0 0 0 0 G 0; 0 0 0 0 0 G];
% 
% S6    = C*E6;
% Cauchy = [S6(1) S6(6) S6(5); S6(6) S6(2) S6(4); S6(5) S6(4) S6(3)];

%This is what I've most recently used... should probably just delete that
%other stuff
% E  = (1/2)*(F'*F-eye(3));
% Cauchy = lambda*trace(E)*eye(3) + 2*mu*E;
    
end