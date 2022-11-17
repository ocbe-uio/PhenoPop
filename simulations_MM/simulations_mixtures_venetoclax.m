function [tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_venetoclax(dose_veneto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function to simulate Venetoclax therapy in three MM patients (MM1420,
%  MM195, MM36) where DECIPHER estimated 2 different clones with 
%  differential responses to the drug. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


init_cell_numb=1.e13; % initial number of cells
Tmax=100;
Tmax2=3000;

%%  Patient MM1420

% Loading results from DECIPHER
results_MM1420_veneto=load('negative_loglikelihood_values-MM1420-drug-1.mat');
num_clones=2; %obtained through the elbow method
param_MM1420_veneto=results_MM1420_veneto.x_finals_temp(:,num_clones,:);

%parameters clone1 obtained by DECIPHER
alpha1=param_MM1420_veneto(1);
b1=param_MM1420_veneto(2);
E1=param_MM1420_veneto(3);
n1=param_MM1420_veneto(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_veneto); % units 1/hour
r1=24*r1_hours; % units 1/day
p1=1.e-13;

%parameters clone2 obtained by DECIPHER
alpha2=param_MM1420_veneto(6);
b2=param_MM1420_veneto(7);
E2=param_MM1420_veneto(8);
n2=param_MM1420_veneto(9);
r2_hours=growth_rate(alpha2,b2,E2,n2,dose_veneto); %units 1/hours
r2=24*r2_hours; %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07; %unit 1/day
 
% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

%initial conditions
prop_clone1=param_MM1420_veneto(5); %proportion of clone1 obtained from DECIPHER
prop_clone2=1-prop_clone1; 
pc_1_init=prop_clone1*init_cell_numb;
pc_2_init=prop_clone2*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

% solving the ode system 
time=[0,100]; 
[tMM1420,x_solMM1420] = ode45(@(t,x) MM_2clone_sys(t,x,teta1),time,x0); 

%%  Patient MM195

% Loading results from DECIPHER
results_MM19520_veneto=load('negative_loglikelihood_values-MM19520-drug-1.mat');
num_clones=2; %obtained through the elbow method
param_MM19520_veneto=results_MM19520_veneto.x_finals_temp(:,num_clones,:);


%parameters clone1 obtained by DECIPHER
alpha1=param_MM19520_veneto(1);
b1=param_MM19520_veneto(2);
E1=param_MM19520_veneto(3);
n1=param_MM19520_veneto(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_veneto); %units 1/hour
r1=24*r1_hours;  %units 1/day
p1=1.e-13;

%parameters clone2 obtained by DECIPHER
alpha2=param_MM19520_veneto(6);
b2=param_MM19520_veneto(7);
E2=param_MM19520_veneto(8);
n2=param_MM19520_veneto(9);
r2_hours=growth_rate(alpha2,b2,E2,n2,dose_veneto); %units 1/hour
r2=24*r2_hours; %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07;

% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

%initial conditions
prop_clone1=param_MM19520_veneto(5); %proportion of clone1 obtained from DECIPHER
prop_clone2=1-prop_clone1;
pc_1_init=prop_clone1*init_cell_numb;
pc_2_init=prop_clone2*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

% solving the ode system 
time=[0,100]; 
[tMM195,x_solMM195] = ode45(@(t,x) MM_2clone_sys(t,x,teta1),time,x0);

%%  patient MM36

% Loading results from DECIPHER
results_MM3620_veneto=load('negative_loglikelihood_values-MM3620-drug-1.mat');
num_clones=2; %obtained through the elbow method
param_MM3620_veneto=results_MM3620_veneto.x_finals_temp(:,num_clones,:);

%parameters clone1
alpha1=param_MM3620_veneto(1);
b1=param_MM3620_veneto(2);
E1=param_MM3620_veneto(3);
n1=param_MM3620_veneto(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_veneto); %unit  1/hour
r1=24*r1_hours; %units 1/day
p1=1.e-13;

%parameters clone2
alpha2=param_MM3620_veneto(6);
b2=param_MM3620_veneto(7);
E2=param_MM3620_veneto(8);
n2=param_MM3620_veneto(9);
r2_hours=growth_rate(alpha2,b2,E2,n2,dose_veneto); %unit 1/hour
r2=24*r2_hours; %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13; 
d3= 0.07; %unit 1/day

% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

%initial conditions
prop_clone1=param_MM3620_veneto(5); %proportion of clone1 obtained from DECIPHER
prop_clone2=1-prop_clone1; 
pc_1_init=prop_clone1*init_cell_numb;
pc_2_init=prop_clone2*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];

% solving the ode system 
time=[0,Tmax2]; 
[tMM36,x_solMM36] = ode45(@(t,x) MM_2clone_sys(t,x,teta1),time,x0); 


end
