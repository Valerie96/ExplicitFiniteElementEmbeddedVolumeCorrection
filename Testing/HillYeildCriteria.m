%Hill Ansisotropic Yeild
mult=1.99;
s11_=1.2*mult; s22_=1.2*mult; s33_=1.2;
s12_=1.8; s13_=1.8; s23_=1.8;

s0=1.2;

R11=s11_/s0; R22=s22_/s0; R33=s33_/s0; 
R12=s12_*sqrt(3)/s0; R13=s13_*sqrt(3)/s0; R23=s23_*sqrt(3)/s0;

F=0.5*(1/R22^2 + 1/R33^2 -1/R11^2);
G=0.5*(1/R33^2 + 1/R11^2 -1/R22^2);
H=0.5*(1/R11^2 + 1/R22^2 -1/R33^2);
L=3/(2*R23^2);
M=3/(2*R13^2);
N=3/(2*R12^2);

s12=0; s13=0; s23=0;

f=@(s11,s22,s33) sqrt(F*(s22-s33).^2+G*(s33-s11).^2+H*(s11-s22).^2+2*L*s23.^2+2*M*s13.^2+2*N*s12.^2);
f2=@(s11,s22,s33) F*(s22-s33).^2+G*(s33-s11).^2+H*(s11-s22).^2+2*L*s23.^2+2*M*s13.^2+2*N*s12.^2;

[s1,s2,s3]=meshgrid([0:0.1:10*s11_]);

figure, hold on, view(3), grid on
isosurface(s1,s2,s3,f2(s1,s2,s3),1)
xlabel('s11');
ylabel('s22');
zlabel('s33');
hold off