function L=ObjectiveFunctionEven(Params,DATA,concvec,timevec,NR)
    %Params are model parameters. DATA is cell counts at each dependent
    %variable condition. ConcVector is list of concentration considered.
    %TimeVector is list of times, and NReps are number of replicates at each
    %condition.
    X1=Params(1);
    Param1=Params(2:5); %alpha,b,E,n pop.1
    Param2=Params(6:9); %alpha,b,E,n pop.2
    sig=Params(10);
    if length(Params) > 10
        error('There are 11 parameters instead of 10 being passed to ObjectiveFunction')
    end
    
    % Param1=[alpha1 b1 E1 n1];
    % Param2=[alpha2 b2 E2 n2];
    NC=length(concvec);
    NT=length(timevec);
    
    sum_resid = 0;
    %X1 = 0.2;
    for i=2:NT%time
        for j=1:NC %concentration
            for k=1:NR % replicate
                
                resid = (DATA(k,j,i) - DATA(k,j,1)*X1*popfunc(Param1,concvec(j),timevec(i))-DATA(k,j,1)*(1-X1)*popfunc(Param2,concvec(j),timevec(i)));
                sum_resid = sum_resid + resid.^2 / (2*sig^2);
            end
        end
    end
    
    
    L = sum_resid + NC*NT*NR*log(sig);
    end
    
    