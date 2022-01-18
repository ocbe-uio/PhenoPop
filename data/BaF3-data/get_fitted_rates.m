function fitted_rates = get_fitted_rates( data , scale)
%this function returns an array of dim( NR, NC ). Each column records the 
%exponential growth rate fit by matlab's fit function with 'exponential'
%method. Matlab fits $\beta e^{\alpha x}$, we fix $\beta$ to be the initial
%time point of the data and fits $\alpha$. 

[NR, NC, NT] = size(data) ;

rep_vec = 1:NR;
conc_vec = scale * [0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6) ];
time_vec = (0:NT-1)*3;


%% Curve Fit

fitted_rates = zeros(NR,NC);
    for conc = 1:NC
        col = rand(1,3);

        for rep = 1 : NR
        cell_count_vec = [];
        for time = 1:NT
            cell_count_vec = [ cell_count_vec;  data(rep,conc,time) ] ;
        end
        
        
        %delete NaN Values and corresponding indicies
        index = isnan(cell_count_vec);
        cell_count_vec(index) = [];
        time = time_vec';
        time(index) = [];
        
        %if entire replicate, conc time series is corrupted, pass
        if isempty(cell_count_vec) == 1
            continue
        end

        y0 = cell_count_vec(1);
        x0 = time(1);
        g = @(p,x)y0*exp(p*(x-x0));
        mdl = fit(time, cell_count_vec, g, 'startpoint', 0.05);
        fitted_rates(rep,conc) = mdl.p;
%         col = rand(1,3);
% 
%         plot(time, cell_count_vec, 'x', 'color', col)
%         hold on
%         plot(time, mdl(time), 'color', col)
% 
    end
end
