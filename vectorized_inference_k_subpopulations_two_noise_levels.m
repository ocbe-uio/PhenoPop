function [fval, x_final]=vectorized_inference_k_subpopulations_two_noise_levels(DATA,no_populations,LML,Conc,Time)

num_optim=500;
NR = 4;
ConcT = 1; %thresholds for sigma
TimeT = 24;
%Conc = [0  0.1 1 10];
%Time = [0 24 48 72];
NP=no_populations;

alpha_ub = 0.10;
b_ub = 1;
E_ub = 50;
n_ub = 50;
if NP>2
    mix_ub=1;
else
    mix_ub=0.5;
end
sig_ub = 5500;


lb= zeros(5*NP + 1,1);
ub=sig_ub*ones(size(lb));

ub(1:5:5*NP-4)=alpha_ub; % alpha
ub(2:5:5*NP-3)=b_ub;  % b
ub(3:5:5*NP-2)=E_ub; % E
ub(4:5:5*NP-1)=n_ub; % n
ub(5:5:5*NP-5)=mix_ub; % Mixture parameter
%lb(3)=0.9;lb(7)=0.9;
lb(2)=0.7;
lb(4)=0.0025;
if NP==2
lb(9)=0.0025;
lb(7)=0.7;
end
% Ensure sum inferred mixture to be < 1
A_inequality=zeros(1,length(ub));
A_inequality(5:5:5*NP-5) = 1;
b_inequality = 1 - LML;


options = optimoptions(@fmincon,'MaxFunctionEvaluations',6990,'MaxIterations',6990,'Display','off');

f=@(x)vectorized_objective_function_k_subpop_and_two_noise_levels(x,no_populations,DATA,Conc,Time, NR, ConcT, TimeT);
fval=inf;
clear xx;
clear ff;
clear x_final;
%rng(42); %set seed to start iterations the same place for every version of R
for nn=1:num_optim
%     if mod(nn,10)==0
%         nn
%     end
    x0=rand(length(ub),1).*(ub-lb) + lb;
    %x0(3:5:5*no_populations-2)=10.^(x0(3:5:5*no_populations-2));    %use when logsampling E
    %lb(3:5:5*no_populations-2)=1e-6; %To avoid Inf or NaN in log() terms we require positive E values
    if no_populations==0
        [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    else
        [xx,ff]=fmincon(f,x0,A_inequality,b_inequality,[],[],lb,ub,[],options);
    end
    if ff<fval
        x_final=xx;
        fval=ff;
    end
end

%       0.026883
%       0.83543
%        822.98
%        0.5145
%       0.19083
%      0.023541
%       0.46274
%        111.29
%        1.5622
%        1040.8
%        375.14
% % Count the number of High and Low noise levels
% Hi_Lo_indicators = ones(NR,NC,NT-1);
% N_hi = sum(Hi_Lo_indicators(:,Conc <= ConcT,Time(2:NT) >= TimeT), 'all');
% N_low = NR*NC*(NT-1) - N_hi;