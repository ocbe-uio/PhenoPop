function plot_linear_timecourse(data,time,name)
%returns plots showing cell counts against time on a log-y scale. Data
% must be either size(14, NC, NT) or size(7, NC, NT). Name is title of the
% plot. 

[NR, NC, NT]=size(data);
if NR == 14
    for conc = 1:NC
        figure
        for rep = 1:NR/2        %make two plots to easily distinguish replicates
            semilogy( time, squeeze(data(rep, conc,:) ), 'linewidth',2 )
            hold on
        end
        title([ name, ' conc ', num2str(conc)]);
        legend('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'location', 'best') 
        hold off
        
        figure
        for rep = (NR/2 + 1):NR
            semilogy( time, squeeze(data(rep, conc,:) ), 'linewidth',2 )
            hold on
        end
        title([ name, ' conc ', num2str(conc)]);
        legend('rep8', 'rep9', 'rep10', 'rep11', 'rep12', 'rep13', 'rep14', 'location', 'best')        
    end
    
elseif NR==7
    for conc = 1:NC
        figure
        for rep = 1:NR
            semilogy( time, squeeze(data(rep, conc,:) ), 'linewidth',2 )
            hold on
            %        ylim([0, 12000])
        end
        title([ name, ' conc ', num2str(conc)]);
        legend('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'location', 'best') 
        hold off
    end
end