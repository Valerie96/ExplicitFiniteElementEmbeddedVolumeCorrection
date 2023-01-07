%% Data vs Flagshyp and Abaqus
file1="ACISpeed2RussellTensile-Half_5000Fibers7";
name1f = "FlagshypACI";
FLAG_1 = ReadFlagshypOutputFile(file1,'jf'); 

suffix = ' ';

[AbqOneHost, AbqETruss, AbqEOne]  = ReadAbaqus_excel(strcat('Abaqus_xlsx/RussellTensile-Half_5000Fibers7_discritized',suffix));
%Russell Tensile 1-5 Nodes: 53, 640 Elements: 181
Abq53=23; Abq640=100; Abq181=1;

C=readcell(strcat("RussellTensionExperiment-10-2", '.xltx'));
Data= cell2mat(C(2:end,:));
Ex_Strain=Data(:,1); Ex_Stress=Data(:,2);
graphsize=[100 100 800 400];
name1 = "Experimental Data";
name1a = "Abaqus Embedded Element";

boundnodes=[631	632	633	634	635	636	637	638	639	640	641	642	643	644	645	646	647	648	649	650	651	652	653	654	655	656	657	658	659	660];
Abn=find(AbqOneHost.nodes==boundnodes(1));
AbRF = AbqOneHost.RF(:,1,Abn);
FRF = FLAG_1.RF(:,1,boundnodes(1));
for i=2:length(boundnodes)
    Abn=find(AbqOneHost.nodes==boundnodes(i));
    AbRF = AbRF+AbqOneHost.RF(:,1,Abn);
    FRF = FRF+FLAG_1.RF(:,1,boundnodes(i));
end

l0=25E-3;
A0=(4E-3)*(6E-3);
%%
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(Ex_Strain,Ex_Stress,'k','DisplayName',name1,'LineWidth',4)
plot(AbqOneHost.U(:,1,Abq640)/l0,AbRF*10^-6/A0,'bo','DisplayName',name1a);
plot(FLAG_1.Disp(:,1,640)/l0,FRF*10^-6/A0,'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Tension Stress vs Strain",suffix));
xlabel("Strain (m/m)");
ylabel("Stress (MPa)");
xlim([0 0.005]);
legend('show');


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.U(:,1,Abq640),'bo','DisplayName',name1a);
plot(FLAG_1.time,FLAG_1.Disp(:,1,640),'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Tension Stress vs Strain",suffix));
xlabel("Time (s)");
ylabel("Disp");
legend('show');
%%
PlotEnergy([AbqEOne.time, AbqEOne.KE],[FLAG_1.Etime, FLAG_1.KE], name1a, name1f, 'Tension - Kinetic Energy')
PlotEnergy([AbqEOne.time, AbqEOne.IE],[FLAG_1.Etime, FLAG_1.IE], name1a, name1f, 'Tension - Internal Energy')
PlotEnergy([AbqEOne.time, -AbqEOne.WK],[FLAG_1.Etime, FLAG_1.WK], name1a, name1f, 'Tension - External Work')
PlotEnergy([AbqEOne.time, AbqEOne.ETOTAL],[FLAG_1.Etime, FLAG_1.ET], name1a, name1f, 'Tension - Total Energy')

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