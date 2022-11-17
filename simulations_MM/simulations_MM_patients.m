%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Script to simulate treatments with Dexamethasone, Selexinor, Meflufen
%  and Venetoclax in three MM patient (MM1420,MM195,MM36) using the
%  clones and corresponding drug sensitivity estimated by DECIPHER.
%
%  Author: Alvaro Köhn-Luque   Date: Nov 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% Dexamethasone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulations
dose_dexa=5; %simulated dexamethasone dose used for all three patients
[tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_dexamethasone(dose_dexa);
Tmax=100;
Tmax2=1000;

%% Plot simulations
fig_dexa=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig_dexa,'defaultAxesColorOrder',[left_color; right_color]);

%Plot patient MM1420
subplot(3,1,1)
yyaxis left 
plot(tMM1420,x_solMM1420(:,3),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM1420,x_solMM1420(:,2),'Color','#813374','LineWidth',2);
hold on
plot(tMM1420,x_solMM1420(:,1),'Color','#c772b9','LineWidth',2);
axis([0 Tmax 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

%legend( 'Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM1420')

%Plot patient MM195
subplot(3,1,2)
yyaxis left
plot(tMM195,x_solMM195(:,2),'k','LineWidth',2);
axis([0 Tmax2 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM195,x_solMM195(:,2),'Color','#72c7b9','LineWidth',2);
hold on
plot(tMM195,x_solMM195(:,1),'Color','#338174','LineWidth',2);
axis([0 Tmax2 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM195')

%Patient MM36
subplot(3,1,3)
yyaxis left
plot(tMM36,x_solMM36(:,2),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM36,x_solMM36(:,2),'Color','#c7c757','LineWidth',2);
hold on
plot(tMM36,x_solMM36(:,1),'Color','#737326','LineWidth',2);
axis([0 Tmax 0 1.1e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM36')

set(gcf,'color','w')
sgtitle('Dexamethasone') ;


%%%%%%%%%%%%%%%%%%%%%%%% Selinexor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulations
dose_seli=0.5; %simulated selinexor dose used for all three patients
[tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_selinexor(dose_seli);
Tmax=100;

%% Plot simulations
fig_seli=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig_seli,'defaultAxesColorOrder',[left_color; right_color]);

%Plot patient MM1420
subplot(3,1,1)
yyaxis left 
plot(tMM1420,x_solMM1420(:,3),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM1420,x_solMM1420(:,2),'Color','#813374','LineWidth',2);
hold on
plot(tMM1420,x_solMM1420(:,1),'Color','#c772b9','LineWidth',2);
axis([0 Tmax 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

%legend( 'Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM1420')

%Plot patient MM195
subplot(3,1,2)
yyaxis left
plot(tMM195,x_solMM195(:,2),'k','LineWidth',2);
axis([0 Tmax 0 4])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM195,x_solMM195(:,2),'Color','#72c7b9','LineWidth',2);
hold on
plot(tMM195,x_solMM195(:,1),'Color','#338174','LineWidth',2);
axis([0 Tmax 0 5.5e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM195')

%Patient MM36
subplot(3,1,3)
yyaxis left
plot(tMM36,x_solMM36(:,2),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM36,x_solMM36(:,2),'Color','#c7c757','LineWidth',2);
hold on
plot(tMM36,x_solMM36(:,1),'Color','#737326','LineWidth',2);
axis([0 Tmax 0 1.1e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM36')

set(gcf,'color','w');
sgtitle('Selinexor')

%%%%%%%%%%%%%%%%%%%%%%%% Meflufen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations
dose_meflu=0.5; %simulated dexamethasone dose used for all three patients
[tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_meflufen(dose_dexa);
Tmax=100;
Tmax2=100;
%% Plot simulations
fig_meflu=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig_meflu,'defaultAxesColorOrder',[left_color; right_color]);

%Plot patient MM1420
subplot(3,1,1)
yyaxis left 
plot(tMM1420,x_solMM1420(:,3),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM1420,x_solMM1420(:,2),'Color','#813374','LineWidth',2);
hold on
plot(tMM1420,x_solMM1420(:,1),'Color','#c772b9','LineWidth',2);
axis([0 Tmax 0 1.1e13])
xlabel('time (days)')
ylabel('cell number')

%legend( 'Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM1420')

%Plot patient MM195
subplot(3,1,2)
yyaxis left
plot(tMM195,x_solMM195(:,2),'k','LineWidth',2);
axis([0 Tmax2 0 4])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM195,x_solMM195(:,2),'Color','#72c7b9','LineWidth',2);
hold on
plot(tMM195,x_solMM195(:,1),'Color','#338174','LineWidth',2);
axis([0 Tmax2 0 5.5e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM195')

%Patient MM36
subplot(3,1,3)
yyaxis left
plot(tMM36,x_solMM36(:,2),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
%plot(tMM36,x_solMM36(:,2),'Color','#c7c757','LineWidth',2);
hold on
plot(tMM36,x_solMM36(:,1),'Color','#737326','LineWidth',2);
axis([0 Tmax 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

%legend('Total M-protein','MM cells clone 1')
title('patient MM36')

set(gcf,'color','w');
sgtitle('Meflufen')

%%%%%%%%%%%%%%%%%%%%%%%% Venetoclax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations
dose_veneto=2; %simulated venetoclax dose used for all three patients
[tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_venetoclax(dose_veneto);
Tmax=100;
Tmax2=3000;

%% Plot simulations
fig_veneto=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig_veneto,'defaultAxesColorOrder',[left_color; right_color]);

%Plot patient MM1420
subplot(3,1,1)
yyaxis left 
plot(tMM1420,x_solMM1420(:,3),'k','LineWidth',2);
axis([0 Tmax 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM1420,x_solMM1420(:,1),'Color','#813374','LineWidth',2);
hold on
plot(tMM1420,x_solMM1420(:,2),'Color','#c772b9','LineWidth',2);
axis([0 Tmax 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

legend( 'Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM1420')

%Plot patient MM195
subplot(3,1,2)
yyaxis left
plot(tMM195,x_solMM195(:,3),'k','LineWidth',2);
axis([0 Tmax 0 5.5])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM195,x_solMM195(:,1),'Color','#338174','LineWidth',2);
hold on
plot(tMM195,x_solMM195(:,2),'Color','#72c7b9','LineWidth',2);
axis([0 Tmax 0 5.e13])
xlabel('time (days)')
ylabel('cell number')

legend('Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM195')

%Patient MM36
subplot(3,1,3)
yyaxis left
plot(tMM36,x_solMM36(:,3),'k','LineWidth',2);
axis([0 Tmax2 0 1.1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM36,x_solMM36(:,2),'Color','#737326','LineWidth',2);
hold on
plot(tMM36,x_solMM36(:,1),'Color','#c7c757','LineWidth',2);
axis([0 Tmax2 0 1.2e13])
xlabel('time (days)')
ylabel('cell number')

legend('Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM36')

set(gcf,'color','w');
sgtitle('Venetoclax')
