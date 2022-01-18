
function pop = vectorized_popfunc(Param,Conc,Time)
    % returns a (1,NC,NT) array 
    %%  Gets population size from model using Hill growth rate function parametrization
    %Param = parameters vector
    %Param(1) alpha
    %Param(2) b
    %Param(3) E
    %Param(4) n
    
    %Conc = drug concentration vector
    %Time = time vector
    
    rate = Param(1) + log( Param(2) + (1-Param(2)) ./ (1+(Conc./Param(3)).^Param(4)) ); %Vector of the same size as Conc
    rep_rate = repmat(rate,1,1,length(Time)); % (1,NC,NT) array
    twist_time = repmat(Time,1,1,length(rate));
    rep_time = permute(twist_time, [1 3 2]);
    pop = exp(rep_time .* rep_rate); % (1,NC,NT) array
    
    end