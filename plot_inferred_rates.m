% Plot rates
%x=linspace(0,5,1000);
x = zeros(1,1000);
x(2:1000) = logspace(-8,-4.86,999);
% Experiment concentrations tested: 
Conc = [0, 31.25*10^-9, 62.5*10^-9, 125*10^-9, 250*10^-9, 375*10^-9, 500*10^-9, 1.25*10^-6, 2.5*10^-6, 3.75*10^-6, 5*10^-6];
N_c = length(Conc);
% Load growth rates found by fitting a slope after taking logarithm of the data
load('./data/dagim_monoclonal/monoclonal_fitted_rates.mat')
fitted_rates_sens500bf
fitted_rates_sens1000bf
fitted_rates_res250bf
fitted_rates_res500bf
fitted_rates_sens500gfp
fitted_rates_sens1000gfp
% Load x finals for inferred rates
load('./data/dagim_monoclonal/monoclonal_x_finals.mat')
%x_final_resistant_250_bf = [0.0461777799165153, 0.000312809399063200, 5.36205275472425e-05, 1.70386099527774]';
%x_final_resistant_500_bf = [0.0439695625181606, 0.00238797357639136, 2.38008734231105e-05, 2.75929869263459]';
%x_final_sensitive_500_bf = [0.0363514973260081, 0.955546663771906, 2.92613885287474e-07, 1.28750067914371]';
%x_final_sensitive_1000_bf = [0.0326113881043157, 0.956838077216945, 4.76581418632678e-07, 1.53401803990261]';
%x_final_sensitive_500_gfp = [0.0370776357879360, 0.852396667612131, 1.22745359434029e-06, 1.17814849337865]';
%x_final_sensitive_1000_gfp = [0.0374774900439094, 0.510891541383781, 1.05801920680056e-05, 0.923224082416028]';

inferred_rate_sens_500_bf = ratefunc(x_final_sensitive_500_bf,x);
inferred_rate_sens_1000_bf = ratefunc(x_final_sensitive_1000_bf,x);
inferred_rate_res_250_bf = ratefunc(x_final_resistant_250_bf,x);
inferred_rate_res_500_bf = ratefunc(x_final_resistant_500_bf,x);
inferred_rate_sens_500_gfp = ratefunc(x_final_sensitive_500_gfp,x);
inferred_rate_sens_1000_gfp = ratefunc(x_final_sensitive_1000_gfp,x);

newcolors0 = [
    0 0 0.726562
    0 0 0.726562
    0.477504 0.821444 0.318195
    0.477504 0.821444 0.318195
    0.792968 0 0
    0.792968 0 0
];

figure
colororder(newcolors0)
semilogx(x+10^(-8),inferred_rate_sens_500_bf, '--', 'linewidth', 2) %, '-o')
hold on 
semilogx(x+10^(-8),inferred_rate_sens_1000_bf, '-', 'linewidth', 2) %, '-+')
semilogx(x+10^(-8),inferred_rate_sens_500_gfp, '--', 'linewidth', 2) %, '-^')
semilogx(x+10^(-8),inferred_rate_sens_1000_gfp, '-', 'linewidth', 2) %, '-p')
semilogx(x+10^(-8),inferred_rate_res_250_bf, '--', 'linewidth', 2) %, '->')
semilogx(x+10^(-8),inferred_rate_res_500_bf, '-', 'linewidth', 2) %, '-v')
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8.3) 10^(-5.1)])
%xlim([1.01 1.010005])
for i=2:length(Conc)
    line([Conc(i)+10^(-8), Conc(i)+10^(-8)], [-0.25, 0.1], 'Color', 'k');
end

semilogx(Conc+10^(-8), fitted_rates_sens500bf, 'xk-')
semilogx(Conc+10^(-8), fitted_rates_sens1000bf, 'xk-')
semilogx(Conc+10^(-8), fitted_rates_res250bf, 'xk-')
semilogx(Conc+10^(-8), fitted_rates_res500bf, 'xk-')
semilogx(Conc+10^(-8), fitted_rates_sens500gfp, 'xk-')
semilogx(Conc+10^(-8), fitted_rates_sens1000gfp, 'xk-')

%xticks(Conc) %logspace(-8,-6,3))
%xticklabels(num2str(Conc))
%xticklabels({['0', '10^-7', '10^-6']})
%yticks([-1 -0.8 -0.2 0 0.2 0.8 1])

title("Inferred growth rates")
legend('Sensitive 500 BF', 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'Observed concentrations', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates.png")

%figure
%semilogx(Conc(2:N_c), fitted_rates_sens500bf(2:N_c), 'x-')
%hold on
%semilogx(Conc(2:N_c), fitted_rates_sens1000bf(2:N_c), 'x-')
%semilogx(Conc(2:N_c), fitted_rates_res250bf(2:N_c), 'x-')
%semilogx(Conc(2:N_c), fitted_rates_res500bf(2:N_c), 'x-')
%semilogx(Conc(2:N_c), fitted_rates_sens500gfp(2:N_c), 'x-')
%semilogx(Conc(2:N_c), fitted_rates_sens1000gfp(2:N_c), 'x-')


vectorized_popfunc(x_final_sensitive_500_bf, Conc, 1);
ratefunc(x_final_sensitive_500_bf, Conc);
