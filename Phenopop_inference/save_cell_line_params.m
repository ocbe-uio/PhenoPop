function out = save_cell_line_params()
%              alpha   b      E    n
%% % high_dE = [0.0275 0.983  0.01 9.1344]; % dE = 1.29
%% high_dE =   [0.0275 0.983  0.7  9.1344]; % dE = 0.6
%% medium_dE = [0.0275 0.983  1.0  9.1344]; % dE = 0.3
%% low_dE =    [0.0275 0.983  1.2  9.1344]; % dE = 0.1
%% fixed =     [0.0275 0.983  1.3  9.1344]; %Resistant

load MIN_BETA.mat
resistant_cell_line = min_beta_E_1975         %[0.0275 0.9877 1.3144  9.1344]; % resistant
original_sensitive_cell_line = min_beta_E_827 %[0.0164 0.9668 10.0000 0.0997]; % sensitive
high_dE = original_sensitive_cell_line

medium_dE = [0.0275 0.9877  0.6  9.1344];
medium_2_dE = [0.0275 0.9877  0.3  9.1344];
low_dE =    [0.0275 0.9877  0.9  9.1344];

% Old setup with med2: 
% Res:  [0.0275 0.9877 1.3144  9.1344]
% Sens: [0.0275 0.9877  0.3  9.1344]

resistant_dagim = [0.045, 0.9, 1.3, 9.1344];
resistant_dagim_5 = [0.035, 0.9, 1.3, 9.1344];
sensitive_dagim = [0.035, 0.9, 0.3, 9.1344];

save CELL_LINE_PARAMETERS.mat resistant_cell_line high_dE medium_dE medium_2_dE low_dE resistant_dagim sensitive_dagim resistant_dagim_5
for i=0:863
   save(strcat('CELL_LINE_PARAMETERS_', int2str(i), '.mat'), 'resistant_cell_line', 'high_dE', 'medium_dE', 'medium_2_dE', 'low_dE', 'resistant_dagim', 'sensitive_dagim', 'resistant_dagim_5')
end
