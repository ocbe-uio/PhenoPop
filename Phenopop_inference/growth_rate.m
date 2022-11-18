function gr = growth_rate(alpha,b,E,n,dose)
% Calculation of the growth rate under treatment

gr=alpha + log(b+(1-b)/(1+(dose/E)^n));

end