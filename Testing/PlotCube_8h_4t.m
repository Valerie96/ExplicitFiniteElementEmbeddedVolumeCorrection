%% Abaqus vs Embedded Abaqus vs   Flagshyp 1 vs Flagshyp 2 vs Flagshyp 3 Energy
cd 'C:\Users\Valerie\Documents\GitHub\ExplicitFiniteElementEmbeddedVolumeCorrection\Testing';
suffix = '';

[AbqOneHost, AbqETruss, AbqEOne]  = ReadAbaqus_excel(strcat('Abaqus_xlsx/FlagshypCube_8h_0t',suffix));
[AbqEHost, AbqETruss, AbqE]  = ReadAbaqus_excel(strcat('Abaqus_xlsx/FlagshypCube_8h_4t',suffix));

graphsize=[100 100 800 400];
name1a = "Abaqus Solid 8h";
name2a = "Abaqus Embedded";

file1=strcat("Cube_8h_0t",suffix);
name1 = "Flagshyp 8h No Truss";

file2=strcat("Cube_8h_4t",suffix);
name2 = "Flagshyp 8h 4t";

file3=strcat("Cube_8h_4t",suffix,"_correct");
name3 = "Flagshyp 8h 4t corrected";

FLAG_1 = ReadFlagshypOutputFile(file1,'job_folder'); 
FLAG_2 = ReadFlagshypOutputFile(file2,'job_folder');
FLAG_3 = ReadFlagshypOutputFile(file3,'job_folder');
% FLAG_2 = ReadFlagshypOutputFile(file2,'job_folder/DataFromPreviousVersion');
% FLAG_3 = ReadFlagshypOutputFile(file3,'job_folder/DataFromPreviousVersion');

PlotEnergy5([AbqEOne.time, AbqEOne.KE],[AbqE.time, AbqE.KE], [FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE],[FLAG_3.Etime, FLAG_3.KE], name1a, name2a,name1, name2,name3,'Kinetic Energy')
PlotEnergy5([AbqEOne.time, AbqEOne.IE],[AbqE.time, AbqE.IE], [FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE],[FLAG_3.Etime, FLAG_3.IE], name1a, name2a,name1, name2,name3,'Internal Energy')
PlotEnergy5([AbqEOne.time, -AbqEOne.WK],[AbqE.time, -AbqE.WK],[FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK],[FLAG_3.Etime, FLAG_3.WK],  name1a, name2a,name1, name2,name3,'External Work')
PlotEnergy5([AbqEOne.time, AbqEOne.ETOTAL],[AbqE.time, AbqE.ETOTAL], [FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET],[FLAG_3.Etime, FLAG_3.ET], name1a, name2a,name1, name2,name3,'Total Energy')

% FLAG_VD3 = ReadFlagshypOutputFileViscDisp(file3);
% PlotEnergy5([AbqEOne.time, AbqEOne.VD], [FLAG_1.Etime, FLAG_1.VD],[FLAG_2.Etime, FLAG_2.VD],[FLAG_VD3.Etime, FLAG_VD3.VD], name1a, name2a,name1, name2,name3,'Viscous Dissipation')
PlotEnergy3([AbqEOne.time, AbqEOne.VD],[AbqE.time, AbqE.VD], [FLAG_2.Etime, FLAG_2.VD], name1a, name2a,name2,'Viscous Dissipation')


IEerror_0t = abs(AbqEOne.IE(end) - FLAG_1.IE(end))/AbqEOne.IE(end)
IEerror_1t = abs(AbqE.IE(end) - FLAG_2.IE(end))/AbqE.IE(end)
IEerror_1tc = abs(AbqEOne.IE(end) - FLAG_3.IE(end))/AbqEOne.IE(end)
%%
F = [1 0.02 0; 0.0 1 0 ;0 0 1-1.912*2E-5];
J     = det(F);  
C     = F'*F;
b     = F*F';  
Ib    = trace(b);     
[V,D] = eig(b); 
nu = 0.3; E=2E11;
K = E/(3*(1-2*nu)); mu = E/(2*(1+nu)); lam = (E*nu/((1+nu)*(1-2*nu)));
C10 = mu/2; D1 = 2/K;
K  = (3*lam+2*mu)/3;
I1  = trace(b);
Cauchy = J^(-5/3)*mu*(b-(1/3)*I1*eye(3)) + K*(J-1)*eye(3);

% E  = (1/2)*(eye(3)-inv(b));
E  = (sqrtm(b)-eye(3));
LogStrain = logm(sqrtm(b));
E=LogStrain;
time = [0:0.01:0.1];

AnalyticDisp3 = ((exp(E(3,3))-1)/2).*ones(1,length(time));
AnalyticDisp2 = ((exp(E(2,2))-1)/2).*ones(1,length(time));
AnalyticDisp1 = (F(1,2)).*ones(1,length(time));

AnalyticStrain1=E(1,1).*ones(1,length(time));
AnalyticStrain2=E(2,2).*ones(1,length(time));

AnalyticStress1=Cauchy(1,1).*ones(1,length(time));
AnalyticStress2=Cauchy(2,2).*ones(1,length(time));
AnalyticStress12=Cauchy(1,2).*ones(1,length(time));

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.RF(:,2,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.RF(:,2,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.RF(:,2,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.RF(:,2,1),'g','DisplayName',name3,'LineWidth',2);
title("External Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
% legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.RF(:,1,27),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.RF(:,1,27),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.RF(:,1,27),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.RF(:,1,27),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.RF(:,1,27),'g','DisplayName',name3,'LineWidth',2);
title("X Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
% legend('show');
%%
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.U(:,1,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.U(:,1,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.Disp(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.Disp(:,1,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.Disp(:,1,1),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticDisp1,'c','DisplayName',"Analytic")
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,1,2),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,1,2),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,1,2),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,1,2),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,1,2),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress1,'c','DisplayName',"Analytic")
title("Host Element 2 XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,1,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,1,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,1,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,1,1),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress1,'c','DisplayName',"Analytic")
title("Host Element 1 XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.RF(:,2,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.RF(:,2,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.RF(:,2,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.RF(:,2,1),'g','DisplayName',name3,'LineWidth',2);
title("External Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,4,2),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,4,2),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,4,2),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,4,2),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,4,2),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress2,'c','DisplayName',"Analytic")
title("Host Element 2 YY Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,4,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,4,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,4,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,4,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,4,1),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress2,'c','DisplayName',"Analytic")
title("Host Element 1 YY Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,2,2),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,2,2),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,2,2),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,2,2),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,2,2),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress12,'c','DisplayName',"Analytic")
title("Host Element 2 XY Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,2,1),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.S(:,2,1),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.HostS(:,2,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,2,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.HostS(:,2,1),'g','DisplayName',name3,'LineWidth',2);
% plot(time,AnalyticStress12,'c','DisplayName',"Analytic")
title("Host Element 1 XY Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');
%%
i=1;
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.U(:,3,i),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.U(:,3,i),'ro','DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.Disp(:,3,i),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.Disp(:,3,i),'r','DisplayName',name1,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.Disp(:,3,i),'g','DisplayName',name1,'LineWidth',2);
title(strcat(int2str(i)," Z Displacement"));
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.RF(:,1,27),'bo','DisplayName',name1a);
plot(AbqE.time,AbqEHost.RF(:,1,27),'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.RF(:,1,27),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.RF(:,1,27),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.RF(:,1,27),'g','DisplayName',name3,'LineWidth',2);
title("X Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
legend('show');
%%
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.S(:,4,1),'bo','DisplayName',name1a);
plot(FLAG_1.time,FLAG_1.HostS(:,4,1),'b','DisplayName',name1,'LineWidth',3);
title("Host YY Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');
%%
for i=1:27
    figure();
    hold on; grid on;
    % fig=gcf; fig.Position=graphsize;
    plot(AbqOneHost.time,AbqOneHost.U(:,2,i),'bo','DisplayName',name1a);
    plot(FLAG_1.time,FLAG_1.Disp(:,2,i),'b','DisplayName',name1,'LineWidth',3);
    title(strcat(int2str(i)," Y Displacement"));
    xlabel("Time (s)");
    ylabel("Displacement (m) ");
    legend('show');
end
%% Function Defs

function PlotEnergy(Data1, Data2, Name1, Name2,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy3(Data1, Data2, Data3, Name1, Name2,Name3,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    plot(Data3(:,1), Data3(:,2),'DisplayName',Name3,'LineWidth',1)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy5(Data1, Data2, Data3, Data4,Data5, Name1, Name2,Name3,Name4,Name5,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'bo','DisplayName',Name1)
    plot(Data2(:,1), Data2(:,2),'ro','DisplayName',Name2)
    plot(Data3(:,1), Data3(:,2),'b','DisplayName',Name3,'LineWidth',3)
    plot(Data4(:,1), Data4(:,2),'r','DisplayName',Name4,'LineWidth',2)
    plot(Data5(:,1), Data5(:,2),'g','DisplayName',Name5,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function [FLAG] = ReadFlagshypOutputFileViscDisp(name)
basedir=strcat('C:/Users/Valerie/Documents/GitHub/FlagshypModified/embeddedelt_edits/job_folder/',name);
f = strcat(basedir,'/VDenergy.dat');
file=fopen(f,'r');
formatSpec = '%e %e';
sizeA = [2 inf ];
TKIE = fscanf(file,formatSpec,sizeA);
fclose(file);

FLAG.Etime =  TKIE(1,:)';
FLAG.VD    =  TKIE(2,:)';
end