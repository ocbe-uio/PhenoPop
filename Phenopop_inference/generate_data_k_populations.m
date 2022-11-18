function [NEWTREAT]=generate_data_k_populations(Conc,Time,NR,MixParams,Params)
%Conc = concentration vector
%Time = time vector
%R = number of replicates
%NoiseL = lower noise level
%NoiseH = higher noise level
%form of MixParams:(p1 p2 p3 .... pk-1) where k is a number of
%subpopulations
%form of Params:(alpha1 b1 E1 n1 .... alpha_k bk Ek nk)
% k_ means (k-1)

NC = length(Conc);
NT = length(Time);
NEWTREAT = zeros(NR, NC, NT); % (ii,jj,kk): Row = replicate, Col = Concentration, Depth = Time

k_=length(MixParams);
MEAN_initial = 1000;
    for ii=1:NR %replicate
        for jj=1:NC %concentration
            for kk=1:NT %time

                for m=1:k_
                    NEWTREAT(ii,jj,kk) = NEWTREAT(ii,jj,kk) + MEAN_initial*MixParams(m)*popfunc(Params(4*m-3:4*m),Conc(jj),Time(kk));
                end
                
                NEWTREAT(ii,jj,kk) = NEWTREAT(ii,jj,kk) + MEAN_initial*(1-sum(MixParams))*popfunc(Params(4*k_+1:4*k_+4),Conc(jj),Time(kk));

                % NEWTREAT(ii,jj,kk) = max(NEWTREAT(ii,jj,kk) + (randn*Noise), 0);

            end

        end
    end

end
