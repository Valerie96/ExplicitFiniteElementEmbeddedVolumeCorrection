%%  Volume Redundancy Guidelines: Dyneema Cube Flagshyp vs Abaqus
cd 'C:/Users/Valerie/OneDrive - The Pennsylvania State University/Research/Projects/GitHub/ExplicitFiniteElementEmbeddedVolumeCorrection/Testing';

file1="VolRedGuidelines_DyneemaCube_0fibers";
name1f = "DyneemaCube noFibers";
FLAG_1 = ReadFlagshypOutputFile(file1,'jf'); 

FLAG_2=FLAG_1; name2f=name1f;
FLAG_3=FLAG_1; name3f=name1f;
% file1="VolRedGuidelines_DyneemaCube_10000fibers";
% name2f = "DyneemaCube_10000fibers";
% FLAG_2 = ReadFlagshypOutputFile(file1,'jf'); 

% file1="ACISpeed3AttwoodCompression-1_1000Fibers7_discritized";
% name3f = "FlagshypACI-SpeedUpdate3";
% FLAG_3 = ReadFlagshypOutputFile(file1,'jf'); 

suffix = '';

[AbqHost1, AbqETruss1, AbqE1]  = ReadAbaqus_excel(strcat("Abaqus_xlsx/DyneemaCube1",suffix));
[AbqHost2, AbqETruss2, AbqE2]  = ReadAbaqus_excel(strcat("Abaqus_xlsx/DyneemaCube1",suffix));

graphsize=[100 100 800 400];
name1a = "Abaqus No Fibers";
name2a = "Abaqus With Fibers";
%%
PlotEnergy5([AbqE1.time, AbqE1.KE],[AbqE2.time, AbqE2.KE], [FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE],[FLAG_3.Etime, FLAG_3.KE], name1a, name2a, name1f, name2f,name3f,'Kinetic Energy')
PlotEnergy5([AbqE1.time, AbqE1.IE],[AbqE2.time, AbqE2.IE], [FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE],[FLAG_3.Etime, FLAG_3.IE], name1a, name2a,name1f, name2f,name3f,'Internal Energy')
PlotEnergy5([AbqE1.time, -AbqE1.WK],[AbqE2.time, -AbqE2.WK],[FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK],[FLAG_3.Etime, FLAG_3.WK],  name1a, name2a,name1f, name2f,name3f,'External Work')
PlotEnergy5([AbqE1.time, AbqE1.ETOTAL],[AbqE2.time, AbqE2.ETOTAL], [FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET],[FLAG_3.Etime, FLAG_3.ET], name1a, name2a,name1f, name2f,name3f,'Total Energy')

%%
boundnodes=[1,  2,  3,  4, 17, 18, 19, 20, 33, 34, 35, 36, 49, 50, 51, 52];
boundelements=[1,  2,  3, 10, 11, 12, 19, 20, 21];

Abn=find(AbqHost1.nodes==boundnodes(1));
AbRF1 = AbqHost1.RF(:,3,Abn);
AbRF2 = AbqHost2.RF(:,3,Abn);
FRF = FLAG_1.RF(:,3,boundnodes(1));
% FRF2 = FLAG_2.RF(:,3,boundnodes(1));
for i=2:length(boundnodes)
    Abn=find(AbqHost1.nodes==boundnodes(i));
    AbRF1 = AbRF1+AbqHost1.RF(:,3,Abn);
    AbRF2 = AbRF2+AbqHost2.RF(:,3,Abn);
    FRF = FRF+FLAG_1.RF(:,3,boundnodes(i));
%     FRF2 = FRF2+FLAG_2.RF(:,3,boundnodes(i));
end

t0=15;
A0=15*15;

Abq11=find(AbqHost1.elements==11);
%%
figure();
hold on; grid on;
plot(-AbqHost1.U(:,3,1)/t0,-AbRF1/A0,'bo','DisplayName',name1a);
plot(-FLAG_1.Disp(:,3,1)/t0,-FRF/A0,'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Compressive Stress vs Strain",suffix));
xlabel("Compressive Strain (m/m)");
ylabel("Stress (MPa)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.LE(:,6,Abq11),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.HostLE(:,6,11),'b','DisplayName',name1f,'LineWidth',3);
title("ZZ Log Strain");
xlabel("Time (s)");
ylabel("Log Strain)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.S(:,6,Abq11),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.HostS(:,6,11),'b','DisplayName',name1f,'LineWidth',3);
title("ZZ Stress");
xlabel("Time (s)");
ylabel("Stress (MPa)");
legend('show');
%%
figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.RF(:,3,1),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.RF(:,3,1),'b','DisplayName',name1f,'LineWidth',3);
title("Z Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.U(:,1,1),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.Disp(:,1,1),'b','DisplayName',name1f,'LineWidth',3);
title(strcat("X Displacement",suffix));
xlabel("Time (s)");
ylabel("Disp (m)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.U(:,3,1),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.Disp(:,3,1),'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Z Displacement",suffix));
xlabel("Time (s)");
ylabel("Disp (m)");
legend('show');
%%
PlotEnergy([AbqE1.time, AbqE1.KE],[FLAG_1.Etime, FLAG_1.KE], name1a, name1f, 'Compression - Kinetic Energy')
PlotEnergy([AbqE1.time, AbqE1.IE],[FLAG_1.Etime, FLAG_1.IE], name1a, name1f, 'Compression - Internal Energy')
PlotEnergy([AbqE1.time, -AbqE1.WK],[FLAG_1.Etime, FLAG_1.WK], name1a, name1f, 'Compression - External Work')
PlotEnergy([AbqE1.time, AbqE1.ETOTAL],[FLAG_1.Etime, FLAG_1.ET], name1a, name1f, 'Compression - Total Energy')
PlotEnergy([AbqE1.time, AbqE1.VD],[FLAG_1.Etime, FLAG_1.VD], name1a, name1f, 'Compression - Viscous Disipation')

%%
PlotEnergy3([AbqE1.time, AbqE1.KE],[FLAG_1.Etime, FLAG_1.KE],[FLAG_2.Etime, FLAG_2.KE], name1a, name1f,name2f,'Compression - Kinetic Energy')
PlotEnergy3([AbqE1.time, AbqE1.IE],[FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE], name1a, name1f,name2f, 'Compression - Internal Energy')
PlotEnergy3([AbqE1.time, -AbqE1.WK],[FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK], name1a, name1f,name2f,'Compression - External Work')
PlotEnergy3([AbqE1.time, AbqE1.ETOTAL],[FLAG_1.Etime, FLAG_1.ET],[FLAG_2.Etime, FLAG_2.ET], name1a, name1f,name2f,'Compression - Total Energy')
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

function PlotEnergy4(Data1, Data2, Data3,Data4, Name1, Name2,Name3,Name4,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',3)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',4)
    plot(Data3(:,1), Data3(:,2),'DisplayName',Name3,'LineWidth',3)
    plot(Data4(:,1), Data4(:,2),'r','DisplayName',Name4,'LineWidth',1)
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