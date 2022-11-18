function [NEWTREAT]=generateDataEven(Conc,Time,NR,MixParam,sensitive,resistant)
%Conc = concentration vector
%Time = time vector
%R = number of replicates
%NoiseL = lower noise level
%NoiseH = higher noise level
%Previously, min_beta_E_827' and min_beta_E_1975' were passed to popfunc
%load MIN_BETA.mat

NC = length(Conc);
NT = length(Time);
NEWTREAT = zeros(NR, NC, NT); % (ii,jj,kk): Row = replicate, Col = Concentration, Depth = Time

MEAN_initial = 1000;
    for ii=1:NR %replicate
        for jj=1:NC %concentration
            for kk=1:NT %time

                NEWTREAT(ii,jj,kk) = MEAN_initial*MixParam*popfunc(sensitive,Conc(jj),Time(kk)) + ...
                    MEAN_initial*(1-MixParam)*popfunc(resistant,Conc(jj),Time(kk));
                
                % NEWTREAT(ii,jj,kk) = max(NEWTREAT(ii,jj,kk) + (randn*Noise), 0);

            end

        end
    end

end



% % 
%plot(concvec, squeeze((NEWTREAT(:,:,5))),'or:');
% 
%plot(concvec, squeeze((NEWTREAT(:,:,6))),'ob:');
%plot(concvec, squeeze((NEWTREAT(:,:,7))),'ok:');
%plot(concvec, squeeze((NEWTREAT(:,:,8))),'og:');
% 
% 
