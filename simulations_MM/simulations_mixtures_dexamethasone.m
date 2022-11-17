function [tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_dexamethasone(dose_dexa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function to simulate Dexamethasone therapy in one MM patient (MM1420),
%  where DECIPHER estimated 2 different clones and two MM patients
%  (MM195, MM36) where DECIPHER estimated a single clone.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init_cell_numb=1.e13; % initial number of cells
Tmax=100;
Tmax2=1000;


%%  Patient MM1420

% Loading results from DECIPHER
results_MM1420_dexa=load('negative_loglikelihood_values-MM1420-drug-4.mat');
num_clones=2; %obtained through the elbow method
param_MM1420_dexa=results_MM1420_dexa.x_finals_temp(:,num_clones,:);

% Parameters clone1 obtained by DECIPHER
alpha1=param_MM1420_dexa(1);
b1=param_MM1420_dexa(2);
E1=param_MM1420_dexa(3);
n1=param_MM1420_dexa(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_dexa); % units 1/hour
r1=24*r1_hours; % units 1/day
p1=1.e-13;

% Parameters clone2 obtained by DECIPHER
alpha2=param_MM1420_dexa(6);
b2=param_MM1420_dexa(7);
E2=param_MM1420_dexa(8);
n2=param_MM1420_dexa(9);
r2_hours=growth_rate(alpha2,b2,E2,n2,dose_dexa); %units 1/hours
r2=24*r2_hours; %units 1/day
p2=1.e-13;

% parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07; %unit 1/day
 
% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

% initial conditions
prop_clone1=param_MM1420_dexa(5); %proportion of clone1 obtained from DECIPHER
prop_clone2=1-prop_clone1; 
pc_1_init=prop_clone1*init_cell_numb;
pc_2_init=prop_clone2*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; pc_2_init; M_pro_init];


% solving the ode system 
time=[0,Tmax]; 
[tMM1420,x_solMM1420] = ode45(@(t,x) MM_2clone_sys(t,x,teta1),time,x0); 

%%  Patient MM195

% Loading results from DECIPHER
results_MM195_dexa=load('negative_loglikelihood_values-MM19520-drug-4.mat');
num_clones=1; %obtained through the elbow method
param_MM195_dexa=results_MM195_dexa.x_finals_temp(:,num_clones,:);

%parameters clone1 obtained by DECIPHER
alpha1=param_MM195_dexa(1);
b1=param_MM195_dexa(2);
E1= param_MM195_dexa(3);
n1=param_MM195_dexa(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_dexa); %units 1/hour
r1=24*r1_hours;  %units 1/day
p1=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07;

% vector with all model parameters
teta1=[r1;p1;r3;d3]; 

%initial conditions
pc_1_init=1*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; M_pro_init];

% solving the ode system 
time=[0,Tmax2]; 
[tMM195,x_solMM195] = ode45(@(t,x) MM_1clone_sys(t,x,teta1),time,x0);

%%  patient MM36

% Loading results from DECIPHER
results_MM36_dexa=load('negative_loglikelihood_values-MM3620-drug-4.mat');
num_clones=1; %obtained through the elbow method
param_MM36_dexa=results_MM36_dexa.x_finals_temp(:,num_clones,:);

%parameters clone1
alpha1=param_MM36_dexa(1);
b1=param_MM36_dexa(2);
E1=param_MM36_dexa(3);
n1=param_MM36_dexa(4);

r1_hours=growth_rate(alpha1,b1,E1,n1,dose_dexa); %unit 1/hour
r1=24*r1_hours; %units 1/day
p1= 1.e-13;

%parameters M-protein
r3= 0.07*1.e-13; 
d3= 0.07; %unit 1/day

% vector with all model parameters
teta1=[r1;p1;r3;d3]; 

%initial conditions
pc_1_init=1*init_cell_numb;
M_pro_init=1;
x0=[pc_1_init; M_pro_init];

% solving the ode system 
time=[0,100]; 
[tMM36,x_solMM36] = ode45(@(t,x) MM_1clone_sys(t,x,teta1),time,x0); 

end
