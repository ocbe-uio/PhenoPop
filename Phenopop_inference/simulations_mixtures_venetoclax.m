%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Script to simulate Venetoclax therapy in three MM patients (MM1420,
%  MM195, MM36) where DECIPHER estimated 2 different clones with 
%  differential responses to the drug. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


init_cell_numb=1.e13 % initial number of cells
dose=2; %drug dose used for all patients

Tmax = 100;
Tmax2 = 3000;

%%  Patient MM1420

%initial conditions
pc_1_init=0.42*init_cell_numb;
pc_2_init=0.58*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

%parameters clone1 obtained by DECIPHER
alpha1=-0.014626333054969;
b1=0.800516341034117;
E1=0.474588337509381;
n1=1.274404720962311;

r1_hours=growth_rate(alpha1,b1,E1,n1,dose); % units 1/hour
r1=24*r1_hours % units 1/day
p1=1.e-13;

%parameters clone2 obtained by DECIPHER
alpha2=-0.0061925409;
b2=0.995727888917675;
E2=565.905337451081;
n2=3.065552447860063;

r2_hours=growth_rate(alpha2,b2,E2,n2,dose); %units 1/hours
r2=24*r2_hours %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07; %unit 1/day
 
% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

% solving the ode system 
time=[0,100]; 
[tMM1420,x_solMM1420] = ode45(@(t,x) MM_venetoclax_sys(t,x,teta1),time,x0); 

%%  Patient MM195

%initial conditions
pc_1_init=0.5*init_cell_numb;
pc_2_init=0.5*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

%parameters clone1 obtained by DECIPHER
alpha1=0.002352212306897;
b1=0.908674363160239;
E1=1257.197507207624;
n1=1.297388511331440;

r1_hours=growth_rate(alpha1,b1,E1,n1,dose); %units 1/hour
r1=24*r1_hours  %units 1/day
p1=1.e-13;

%parameters clone2 obtained by DECIPHER
alpha2=0.00330548442;
b2=0.993327491549934;
E2=44.229024269839;
n2=7.302229703119726;

r2_hours=growth_rate(alpha2,b2,E2,n2,dose); %units 1/hour
r2=24*r2_hours %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07;

% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

% solving the ode system 
time=[0,100]; 
[tMM195,x_solMM195] = ode45(@(t,x) MM_venetoclax_sys(t,x,teta1),time,x0);

%%  patient MM36
%initial conditions
pc_1_init=0.23*init_cell_numb;
pc_2_init=0.77*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

%parameters clone1
alpha1=0.000050469350237;
b1=0.982794370259542;
E1=163.781986574522;
n1=1.435175376044731;

r1_hours=growth_rate(alpha1,b1,E1,n1,dose); %unit  1/hour
r1=24*r1_hours %units 1/day
p1= 1.e-13;

%parameters clone2
alpha2=0.0001490435;
b2=0.759865318430992;
E2=98.242804796922;
n2=1.678194166301342;

r2_hours=growth_rate(alpha2,b2,E2,n2,dose); %unit  1/hour
r2=24*r2_hours %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13; 
d3= 0.07; %unit 1/day

% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

% solving the ode system 
time=[0,3000]; 
[tMM36,x_solMM36] = ode45(@(t,x) MM_venetoclax_sys(t,x,teta1),time,x0); 


%% Plot simulations
fig=figure
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

%Plot patient MM1420
subplot(3,1,1)
yyaxis left 
plot(tMM1420,x_solMM1420(:,3),'k','LineWidth',2);
axis([0 Tmax 0 1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM1420,x_solMM1420(:,1),'Color','#813374','LineWidth',2);
hold on
plot(tMM1420,x_solMM1420(:,2),'Color','#c772b9','LineWidth',2);
axis([0 Tmax 0 6e12])
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
plot(tMM195,x_solMM195(:,1),'Color','#72c7b9','LineWidth',2);
hold on
plot(tMM195,x_solMM195(:,2),'Color','#338174','LineWidth',2);
axis([0 Tmax 0 5.e13])
xlabel('time (days)')
ylabel('cell number')

legend('Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM195')

%Patient MM36
subplot(3,1,3)
yyaxis left
plot(tMM36,x_solMM36(:,3),'k','LineWidth',2);
axis([0 Tmax2 0 1])
xlabel('time (days)')
ylabel('M-protein (g/dL)')

yyaxis right
plot(tMM36,x_solMM36(:,1),'Color','#c7c757','LineWidth',2);
hold on
plot(tMM36,x_solMM36(:,2),'Color','#737326','LineWidth',2);
%axis([0 Tmax2 0 1e12])
xlabel('time (days)')
ylabel('cell number')

legend('Total M-protein','MM cells clone 1','MM cells clone 2')
title('patient MM36')

set(gcf,'color','w');


