%--------------------------------------------------------------------------
% Computes the element vector of global internal forces 
%--------------------------------------------------------------------------
function [T_internal,PLAST,geomJn_1,VolRate,Cauchy,epsilon,CauchyTensor] = element_force_truss(...
          matyp,properties,x_local,X_local,PLAST,GEOM,DAMPING,dt)  
    
%--------------------------------------------------------------------------
% Temporary variables.
%--------------------------------------------------------------------------
L       = norm(X_local(:,2) - X_local(:,1));  
dx      = x_local(:,2) - x_local(:,1);        
l       = norm(dx);                           
n       = dx/l;                                                           
lambda  = l/L;                                
epsilon = log(lambda);

%Bulk Viscosity Damping
% Jn_1=GEOM.Jn_1(ielement);
% eps_dot = (J-Jn_1)/dt;
% b1=DAMPING.b1;
% b2=DAMPING.b2;
% 
%     le=calc_element_size(FEM,GEOM,ielement);
%     Cd=sqrt((lambda + 2*mu)/rho);
%     
%     p1 = rho*b1*le*Cd*eps_dot*eye(3);
%     p2 = rho*(b2*le)^2*abs(eps_dot)*min(0,eps_dot)*eye(3);

switch matyp
    case 2
        %Explicit NeoHooke
        % rho   = properties(1);
        E     = properties(2);
        nu    = properties(3);
        area  = properties(4);
        % ty0   = properties(5);  
        % H     = properties(6);
        lam = (E*nu/((1+nu)*(1-2*nu))); mu = E/(2*(1+nu)); 
        
        % ep    = PLAST.ep;    
        % epbar = PLAST.epbar; 

        J       = lambda^(1-2*nu);
        b = [lambda^2 0 0; 0 lambda^(-2*nu) 0; 0 0 lambda^(-2*nu)];
        s = (1/3)*(3*lam+2*mu)*(J-1)*eye(3) + mu*(J^(-5/3))*(b - (1/3)*trace(b))*eye(3);
        CauchyTensor = s;
        Cauchy = CauchyTensor(1,1);

    case 4
        %Ogden
        % rho   = properties(1);
        alpha1     = properties(2);
        mu1        = properties(3);
        D1         = properties(4);
        area       = properties(5);
        K=2/D1;
        mu = mu1*alpha1/2;
        nu = (3*K-2*mu)/(2*(3*K+mu));
        
        J       = lambda^(1-2*nu);
        %b = [lambda^2 0 0; 0 lambda^(-2*nu) 0; 0 0 lambda^(-2*nu)];
        
        lambda1 = lambda/J^(1/3);
        lambda2 = (lambda^-nu)/J^(1/3);
        lambda3 = (lambda^-nu)/J^(1/3);
        PF1     = 2*mu1/(alpha1*J);
        
        sigma1 = (PF1*(lambda1^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1));
        sigma2 = (PF1*(lambda2^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1));
        sigma3 = (PF1*(lambda3^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1));
        CauchyTensor =[sigma1 0 0; 0 sigma2 0; 0 0 sigma3];
        Cauchy = CauchyTensor(1,1);
end


%--------------------------------------------------------------------------
%Abaqus considers trusses incompressible and therefore does not calculate
%any change in cross sectional area that ought to happen for large strain
%problems. This is the code to calculate that change, which will not be
%used while we are trying to mimic Abaqus

% epsilon = (l-L)/L;
% Cauchy = E*epsilon;
%%area might need to get smaller
%  a = area*L*J/l;

%Update previous Jacobian and element strain rate
%     geomJn_1=J;
%     VolRate = eps_dot;
%--------------------------------------------------------------------------
a = area*L/l;
% a = area;
T          = Cauchy*a;                  
Tb         = T*n;                         
T_internal = [-Tb;Tb];    

    
%Update previous Jacobian and element strain rate
    geomJn_1=1;
    VolRate = 1;

end