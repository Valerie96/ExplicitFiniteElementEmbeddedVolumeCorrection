%Ogden Material Model (from Ritika)
alpha1 = 4.5;
mu1 = 2.5e3;
D1=0.9091e-9;
K=2/D1;

X=[1 1 1;1 0 1;1 1 0;1 0 0;0 1 1;0 0 1;0 1 0;0 0 0]';

F = [1, 0, 0; 1, 1, 0; 0, 0, 1];

J = det(F);
B = F*F';
[V,D] = eig(B);

lambda1 = sqrt(D(1,1))/J^(1/3);
lambda2 = sqrt(D(2,2))/J^(1/3);
lambda3 = sqrt(D(3,3))/J^(1/3);
PF1 = 2*mu1/(alpha1*J);
PFd1 = mu1/J;


LE=zeros(3);
for i=1:3
    LE = LE + log(D(i,i))*V(:,i)*V(:,i)';
end
%Abaqus LEii is divided by 2
LE(1,1)=LE(1,1)/2; LE(2,2)=LE(2,2)/2; LE(3,3)=LE(3,3)/2;

sigmaabaqus1 = (PF1*(lambda1^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1))
sigmaabaqus2 = (PF1*(lambda2^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1))
sigmaabaqus3 = (PF1*(lambda3^alpha1 - (lambda1^alpha1+lambda2^alpha1+lambda3^alpha1)/3) + K*(J-1))


x=F*X;
u=x-X;

%Verified for cubes:
%F = [1.5, 0, 0; 0, 1, 0; 0, 0, 1];
% F = [1.5, 0, 0; 0, 0.7, 0; 0, 0, 0.7];
% F = [1, 0, 0; 0, 0.9, 0; 0, 0, 2.2];

%%
%Mooney Rivlin
C10=-100;
C01=1.2e3;
D1=0.9091e-9;
K = 2/D1;

X=[1 1 1;1 0 1;1 1 0;1 0 0;0 1 1;0 0 1;0 1 0;0 0 0]';

F = [1, 0, 0; 1, 1, 0; 0, 0, 1];
%  F = [1, 0, 0; 0, 0.9, 0; 0, 0, 2.2];

F = [2, 0, 0; 0, 1-0.2924, 0; 0, 0, 1-0.2924];
F = [2, 0, 0; 0, 0.7, 0; 0, 0, 0.7];
J = det(F);
B = F*F';
[V,D] = eig(B);
VT = V';
hydro = K*(J-1);
lambda1 = sqrt(D(1,1))/J^(1/3);
lambda2 = sqrt(D(2,2))/J^(1/3);
lambda3 = sqrt(D(3,3))/J^(1/3);

I1_ = lambda1^2+lambda2^2+lambda3^2;
I2_ = lambda1^2*lambda2^2+lambda2^2*lambda3^2+lambda3^2*lambda1^2;

LE=zeros(3);
for i=1:3
    LE = LE + log(D(i,i))*V(:,i)*V(:,i)';
end
%Abaqus LEii is divided by 2
LE(1,1)=LE(1,1)/2; LE(2,2)=LE(2,2)/2; LE(3,3)=LE(3,3)/2;

LE
x=F*X;
u=x-X;

C1=C10; C2=C01;

Cauchy=(2/J)*(J^(-2/3)*(C1+I1_*C2)*B - J^(-4/3)*C2*B*B) + ((2/D1)*(J-1)-2/(3*J)*(C1*I1_+2*C2*I2_))*eye(3)

% Cauchy=(1/J)*(2*J^(-2/3)*(C1+I1_*C2)*B-2*J^(-4/3)*C2*B*B-(2/3)*(C1*I1_+2*C2*I2_)*eye(3))+2/D1*(J-1)

%My derivation
K               = 2/D1;
I1              = trace(B);
I2              = (1/2)*(trace(B)^2-trace(B*B));
Cauchy          = J^(-5/3)*2*C10*(B-(1/3)*I1*eye(3)) + J^(-7/3)*2*C01*(I1*B-B*B-(2/3)*I2*eye(3)) + K*(J-1)*eye(3)

%Verified for cubes:
% F = [2, 0, 0; 0, 0.7, 0; 0, 0, 0.7];
% F = [1, 0, 0; 0, 0.9, 0; 0, 0, 2.2];


% B_=det(B)^(-1/3)*B;
% [V_ D_]=eig(B_);
% LE_=zeros(3);
% for i=1:3
%     LE_ = LE_ + log(D_(i,i))*V_(:,i)*V_(:,i)';
% end


%% 
% Neo Hooke
C10=1100;
D1=0.9091e-9;
K = 2/D1;
mu = 2*C10;

X=[1 1 1;1 0 1;1 1 0;1 0 0;0 1 1;0 0 1;0 1 0;0 0 0]';

F = [1, 0.1, 0; 0.1, 1, 0; 0, 0, 1-8.155e-3];
%  F = [1, 0, 0; 0, 0.9, 0; 0, 0, 2.2];
J = det(F);
B = F*F';
[V,D] = eig(B);


LE=zeros(3);
for i=1:3
    LE = LE + log(D(i,i))*V(:,i)*V(:,i)';
end
%Abaqus LEii is divided by 2
LE(1,1)=LE(1,1)/2; LE(2,2)=LE(2,2)/2; LE(3,3)=LE(3,3)/2;

LE
x=F*X;
u=x-X;

I1              = trace(B);
Cauchy          = J^(-5/3)*2*C10*(B-(1/3)*I1*eye(3)) + K*(J-1)*eye(3)
