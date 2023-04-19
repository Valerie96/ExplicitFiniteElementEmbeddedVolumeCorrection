
cd 'C:/Users/Valerie/OneDrive - The Pennsylvania State University/Research/Projects/GitHub/ExplicitFiniteElementEmbeddedVolumeCorrection/Testing';

file1="MRTension_1h";
name1 = "Flagshyp 1h";

file2="MRTension_1h";
name2 = "Flagshyp 1h Disp Spec";

FLAG_1 = ReadFlagshypOutputFile(file1,'job_folder'); 
FLAG_2 = ReadFlagshypOutputFile(file2,'job_folder'); 

figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.HostS(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,1,1),'r','DisplayName',name2,'LineWidth',3);
% plot(time,AnalyticStress1,'c','DisplayName',"Analytic")
title("Flagshyp Host Element 1 XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
plot(FLAG_1.HostLE(:,1,1),FLAG_1.HostS(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.HostLE(:,1,1),FLAG_2.HostS(:,1,1),'r','DisplayName',name2,'LineWidth',3);
% plot(time,AnalyticStress1,'c','DisplayName',"Analytic")
title("Host Element 1 XX Stress");
xlabel("Strain");
ylabel("Stress (Pa)");
legend('show');
%%
file1="NHTension_1h";
name1 = "Flagshyp 1h";

file2="MRTension_1h";
name2 = "Flagshyp 1h MR";

FLAG_1 = ReadFlagshypOutputFile(file1,'job_folder'); 
FLAG_2 = ReadFlagshypOutputFile(file2,'job_folder'); 



figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.Disp(:,2,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.Disp(:,2,1),'r','DisplayName',name2,'LineWidth',3);
% plot(time,AnalyticDisp1,'c','DisplayName',"Analytic")
title("Y Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.Disp(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.Disp(:,1,1),'r','DisplayName',name2,'LineWidth',3);
% plot(time,AnalyticDisp1,'c','DisplayName',"Analytic")
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');



figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.HostS(:,1,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.HostS(:,1,1),'r','DisplayName',name2,'LineWidth',3);
% plot(time,AnalyticStress1,'c','DisplayName',"Analytic")
title("Host Element 1 XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');


%%
C10=-100;
C01=1200;
D1=0.9091e-9;
mu = 2*(C10+C01)
K=2/D1; lam=K-2*mu/3