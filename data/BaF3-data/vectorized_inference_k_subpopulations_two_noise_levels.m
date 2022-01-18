function [fval, x_final, N_hi, N_low ]=vectorized_inference_k_subpopulations_two_noise_levels(DATA, Conc, ConcT, TimeT, ub, lb, num_optim, no_populations, A_inequality, b_inequality)

[NR, NC, NT]=size(DATA);

Time = [0 : NT-1]*3;
options = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',5990,'Display','off');

f=@(x)vectorized_objective_function_two_noise_levels_k_subpopulations(x,no_populations,DATA,Conc,Time, ConcT, TimeT, NR);
fval=inf;
clear xx;
clear ff;
clear x_final;
%rng(42); %set seed to start iterations the same place for every version of R
for nn=1:num_optim
    
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

% % Count the number of High and Low noise levels
% Hi_Lo_indicators = ones(NR,NC,NT-1);
% N_hi = sum(Hi_Lo_indicators(:,Conc <= ConcT,Time(2:NT) >= TimeT), 'all');
% N_low = NR*NC*(NT-1) - N_hi;