load('./data/NSCLC-and-BaF3-cell-data/TREATS.mat') %Load data.
Conc=[0 0.1 1 10];
Time=[0 24 48 72];
%calculate growth rates for monoclonal populations.
[fval827, x_final827]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_827,1,0,Conc,Time);
[fval1975, x_final1975]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_1975,1,0,Conc,Time);

% %calculate growth rates for mixture populations.
[fvalMIX11, x_finalMIX11]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_11,2,0,Conc,Time);
[fvalMIX41, x_finalMIX41]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_41,2,0,Conc,Time);
%[fvalMIX14, x_finalMIX14]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_14,2,0,Conc,Time);
%[fvalMIX19, x_finalMIX19]=vectorized_inference_k_subpopulations_two_noise_levels(TREAT_E_19,2,0,Conc,Time);

GR827=x_final827(1:4);
GR1975=x_final1975(1:4);

GRMIX11_1=x_finalMIX11(1:4);
GRMIX11_2=x_finalMIX11(6:9);

GRMIX41_1=x_finalMIX41(1:4);
GRMIX41_2=x_finalMIX41(6:9);

%GRMIX14_1=x_finalMIX14(1:4);
%GRMIX14_2=x_finalMIX14(6:9);

%GRMIX19_1=x_finalMIX19(1:4);
%GRMIX19_2=x_finalMIX19(6:9);

z827 = 0.5*(GR827(1)+growth_rate(GR827(1), GR827(2), GR827(3), GR827(4),10));
z1975 = 0.5*(GR1975(1)+growth_rate(GR1975(1), GR1975(2), GR1975(3), GR1975(4),10));
zMIX11_1=0.5*(GRMIX11_1(1)+growth_rate(GRMIX11_1(1), GRMIX11_1(2), GRMIX11_1(3), GRMIX11_1(4),10));
zMIX11_2=0.5*(GRMIX11_2(1)+growth_rate(GRMIX11_2(1), GRMIX11_2(2), GRMIX11_2(3), GRMIX11_2(4),10));
zMIX41_1=0.5*(GRMIX41_1(1)+growth_rate(GRMIX41_1(1), GRMIX41_1(2), GRMIX41_1(3), GRMIX41_1(4),10));
zMIX41_2=0.5*(GRMIX41_2(1)+growth_rate(GRMIX41_2(1), GRMIX41_2(2), GRMIX41_2(3), GRMIX41_2(4),10));
%zMIX14_1=0.5*(GRMIX14_1(1)+growth_rate(GRMIX14_1(1), GRMIX14_1(2), GRMIX14_1(3), GRMIX14_1(4),10));
%zMIX14_2=0.5*(GRMIX14_2(1)+growth_rate(GRMIX14_2(1), GRMIX14_2(2), GRMIX14_2(3), GRMIX14_2(4),10));
%zMIX19_1=0.5*(GRMIX19_1(1)+growth_rate(GRMIX19_1(1), GRMIX19_1(2), GRMIX19_1(3), GRMIX19_1(4),10));
%zMIX19_2=0.5*(GRMIX19_2(1)+growth_rate(GRMIX19_2(1), GRMIX19_2(2), GRMIX19_2(3), GRMIX19_2(4),10));

gc50_827=GR827(3)*((exp(z827-GR827(1))-1)/(GR827(2)-exp(z827-GR827(1))))^(1/GR827(4));
gc50_1975=GR1975(3)*((exp(z1975-GR1975(1))-1)/(GR1975(2)-exp(z1975-GR1975(1))))^(1/GR1975(4));

gc50_MIX11_1=GRMIX11_1(3)*((exp(zMIX11_1-GRMIX11_1(1))-1)/(GRMIX11_1(2)-exp(zMIX11_1-GRMIX11_1(1))))^(1/GRMIX11_1(4));
gc50_MIX11_2=GRMIX11_2(3)*((exp(zMIX11_2-GRMIX11_2(1))-1)/(GRMIX11_2(2)-exp(zMIX11_2-GRMIX11_2(1))))^(1/GRMIX11_2(4));

gc50_MIX41_1=GRMIX41_1(3)*((exp(zMIX41_1-GRMIX41_1(1))-1)/(GRMIX41_1(2)-exp(zMIX41_1-GRMIX41_1(1))))^(1/GRMIX41_1(4));
gc50_MIX41_2=GRMIX41_2(3)*((exp(zMIX41_2-GRMIX41_2(1))-1)/(GRMIX41_2(2)-exp(zMIX41_2-GRMIX41_2(1))))^(1/GRMIX41_2(4));

%gc50_MIX14_1=GRMIX14_1(3)*((exp(zMIX14_1-GRMIX14_1(1))-1)/(GRMIX14_1(2)-exp(zMIX14_1-GRMIX14_1(1))))^(1/GRMIX14_1(4));
%gc50_MIX14_2=GRMIX14_2(3)*((exp(zMIX14_2-GRMIX14_2(1))-1)/(GRMIX14_2(2)-exp(zMIX14_2-GRMIX14_2(1))))^(1/GRMIX14_2(4));

%gc50_MIX19_1=GRMIX19_1(3)*((exp(zMIX19_1-GRMIX19_1(1))-1)/(GRMIX19_1(2)-exp(zMIX19_1-GRMIX19_1(1))))^(1/GRMIX19_1(4));
%gc50_MIX19_2=GRMIX19_2(3)*((exp(zMIX19_2-GRMIX19_2(1))-1)/(GRMIX19_2(2)-exp(zMIX19_2-GRMIX19_2(1))))^(1/GRMIX19_2(4));



