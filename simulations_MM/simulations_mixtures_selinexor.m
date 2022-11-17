function [tMM1420,x_solMM1420,tMM195,x_solMM195,tMM36,x_solMM36]=simulations_mixtures_selinexor(dose_seli)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Script to simulate Selinexor therapy in one MM patient (MM1420),
%  where DECIPHER estimated 2 different clones and two other patients
%  (MM195, MM36) where DECIPHER estimated a single clone.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init_cell_numb=1.e13; % initial number of cells
Tmax=100;
Tmax2=100;


%%  Patient MM1420

% Loading results from DECIPHER
results_MM1420_seli=load('negative_loglikelihood_values-MM1420-drug-3.mat');
num_clones=2; %obtained through the elbow method
param_MM1420_seli=results_MM1420_seli.x_finals_temp(:,num_clones,:);

%parameters clone1 obtained by DECIPHER
alpha1=param_MM1420_seli(1);
b1=param_MM1420_seli(2);
E1=param_MM1420_seli(3);
n1=param_MM1420_seli(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_seli); % units 1/hour
r1=24*r1_hours; % units 1/day
p1=1.e-13;

%parameters clone2 obtained by DECIPHER
alpha2=param_MM1420_seli(6);
b2=param_MM1420_seli(7);
E2=param_MM1420_seli(8);
n2=param_MM1420_seli(9);
r2_hours=growth_rate(alpha2,b2,E2,n2,dose_seli); %units 1/hours
r2=24*r2_hours; %units 1/day
p2=1.e-13;

%parameters M-protein
r3= 0.07*1.e-13;
d3= 0.07; %unit 1/day
 
% vector with all model parameters
teta1=[r1;p1;r2;p2;r3;d3]; 

%initial conditions
prop_clone1=param_MM1420_seli(5);
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
results_MM195_seli=load('negative_loglikelihood_values-MM19520-drug-3.mat');
num_clones=1; %obtained through the elbow method
param_MM195_seli=results_MM195_seli.x_finals_temp(:,num_clones,:);

%parameters clone1 obtained by DECIPHER
%alpha1=0.002707881756636;
%b1=0.993219048240030;
%E1= 9.999261793388487e+03;
%n1=0.050257476348491;
alpha1=param_MM195_seli(1);
b1=param_MM195_seli(2);
E1=param_MM195_seli(3);
n1=param_MM195_seli(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_seli); %units 1/hour
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
results_MM36_seli=load('negative_loglikelihood_values-MM3620-drug-3.mat');
num_clones=1; %obtained through the elbow method
param_MM36_seli=results_MM36_seli.x_finals_temp(:,num_clones,:);

%parameters clone1
%alpha1=-0.002616609844500;
%b1=0.998378901463355;
%E1=4.204600770400758;
%n1=55.648038136731046;
alpha1=param_MM36_seli(1);
b1=param_MM36_seli(2);
E1=param_MM36_seli(3);
n1=param_MM36_seli(4);
r1_hours=growth_rate(alpha1,b1,E1,n1,dose_seli); %unit 1/hour
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