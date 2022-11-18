function [pop] = popfunc(b,x, T)
%%  Gets population size from model using Hill growth rate function parametrization
%b = parameters vector
%b(1) alpha
%b(2) b
%b(3) E
%b(4) n

%x = drug concentration
%T = time


%x concentration
pop = (exp(T* (b(1)+log(b(2)+ (1-b(2))./(1+(x./b(3)).^b(4))))));
 
end

