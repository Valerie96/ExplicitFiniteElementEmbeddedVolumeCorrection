%%  Volume Redundancy Guidelines: Dyneema Cube Flagshyp vs Abaqus
cd 'C:\Users\Valerie\Documents\GitHub\ExplicitFiniteElementEmbeddedVolumeCorrection\Testing';

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

suffix = ' ';

[AbqHost1, AbqETruss1, AbqE1]  = ReadAbaqus_excel(strcat("Abaqus_xlsx/DyneemaCube",suffix));
[AbqHost2, AbqETruss2, AbqE2]  = ReadAbaqus_excel(strcat("Abaqus_xlsx/DyneemaCube_10000fibers",suffix));

graphsize=[100 100 800 400];
name1a = "Abaqus No Fibers";
name2a = "Abaqus With Fibers";
%%
PlotEnergy5([AbqE1.time, AbqE1.KE],[AbqE2.time, AbqE2.KE], [FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE],[FLAG_3.Etime, FLAG_3.KE], name1a, name2a, name1f, name2f,name3f,'Kinetic Energy')
PlotEnergy5([AbqE1.time, AbqE1.IE],[AbqE2.time, AbqE2.IE], [FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE],[FLAG_3.Etime, FLAG_3.IE], name1a, name2a,name1f, name2f,name3f,'Internal Energy')
PlotEnergy5([AbqE1.time, -AbqE1.WK],[AbqE2.time, -AbqE2.WK],[FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK],[FLAG_3.Etime, FLAG_3.WK],  name1a, name2a,name1f, name2f,name3f,'External Work')
PlotEnergy5([AbqE1.time, AbqE1.ETOTAL],[AbqE2.time, AbqE2.ETOTAL], [FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET],[FLAG_3.Etime, FLAG_3.ET], name1a, name2a,name1f, name2f,name3f,'Total Energy')

%%
boundnodes=[1,2,3,4,5,6,7,8,9,10,11,122,123,124,125,126 ...
127,128,129,130,131,132,243,244,245,246,247,248,249,250,251,252 ...
253,364,365,366,367,368,369,370,371,372,373,374,485,486,487,488 ...
489,490,491,492,493,494,495,606,607,608,609,610,611,612,613,614 ...
615,616,727,728,729,730,731,732,733,734,735,736,737,848,849,850 ...
851,852,853,854,855,856,857,858,969,970,971,972,973,974,975,976 ...
977,978,979,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1211 ...
1212,1213,1214,1215,1216,1217,1218,1219,1220,1221];


boundelements=[1,2,3,4,5,6,7,8,9,10,101,102,103,104,105,106 ...
107,108,109,110,201,202,203,204,205,206,207,208,209,210,301,302 ...
303,304,305,306,307,308,309,310,401,402,403,404,405,406,407,408 ...
409,410,501,502,503,504,505,506,507,508,509,510,601,602,603,604 ...
605,606,607,608,609,610,701,702,703,704,705,706,707,708,709,710 ...
801,802,803,804,805,806,807,808,809,810,901,902,903,904,905,906 ...
907,908,909,910];

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

%Average center stress
%Center elements: 405,505,406,506 Ceter node: 611
Abq611=116;
Abq405=find(boundelements==405);Abq505=find(boundelements==505);
Abq406=find(boundelements==406);Abq506=find(boundelements==506);
Ab1_S33=mean(AbqHost1.S(:,3,Abq405),AbqHost1.S(:,3,Abq505),AbqHost1.S(:,3,Abq406),AbqHost1.S(:,3,Abq506));
%%
figure();
hold on; grid on;
plot(-AbqHost1.U(:,3,Abq611)/t0,-AbRF1/A0,'bo','DisplayName',name1a);
plot(-FLAG_1.Disp(:,3,611)/t0,-FRF/A0,'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Compressive Stress vs Strain",suffix));
xlabel("Compressive Strain (m/m)");
ylabel("Stress (MPa)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.S(:,3,Abq405),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.HostS(:,3,405),'b','DisplayName',name1f,'LineWidth',3);
title("Stress");
xlabel("Time (s)");
ylabel("Stress (MPa)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.RF(:,3,Abq611),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.RF(:,3,611),'b','DisplayName',name1f,'LineWidth',3);
title("X Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
legend('show');

figure();
hold on; grid on;
plot(AbqHost1.time,-AbqHost1.U(:,3,Abq611),'bo','DisplayName',name1a);
plot(FLAG_1.time,-FLAG_1.Disp(:,3,611),'b','DisplayName',name1f,'LineWidth',3);
title(strcat("Displacement",suffix));
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