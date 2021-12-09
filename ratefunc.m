function [rate] = ratefunc(b,x) %alpha b E N

rate = (b(1)+log(b(2)+ (1-b(2))./(1+(x./b(3)).^b(4))));
end

%xlinconc = linspace(0.01,10,100);
% hold on
 % for j=1:4 %replicates
%        semilogx(concvec_fake,expr(:,j),'bx')
%    end
%legend( 'Exponential growth rates from data','Model fit')