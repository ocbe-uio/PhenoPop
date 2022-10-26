% Turn off figures popping up:
set(0,'DefaultFigureVisible','off')

NR = 4;
Conc = [0.0, 0.000005, 0.000010772173450, 0.000023207944168, 0.00005, 0.000107721734502, 0.000232079441681, 0.0005, 0.001077217345016, 0.002320794416806, 0.005, 0.010772173450159, 0.023207944168064, 0.05, 0.107721734501594, 0.232079441680639, 0.5];
Time = [0 12 24 36 48 60 72 84 96]; % hours

% The example data is the simulated "DATA-2-case-1-Noise-50-mix-0.5.csv"
datafile = "./example_data.csv";
[x_finals_temp, f_vals_temp, negative_loglikelihood_values, AIC_values, BIC_values] = PhenoPop("Example data", NR, Conc, Time, datafile, 3, 10); %, 4, 20, [-0.1, 0, 1e-6,   0, 1e-6], [ 0.1, 1,  2e3, 100, 5e4], true, 10, 48);
