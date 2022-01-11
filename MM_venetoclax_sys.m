function dxdt = MM_2drugs_sys(t,x,teta)
%ode system for MM dynamics with 2 clones

%notation of model parameters   
r1=teta(1);
p1=teta(2);
r2=teta(3);
p2=teta(4);
r3=teta(5);
d3=teta(6);


dxdt=zeros(3,1);

dxdt(1) = r1*x(1)*(1/(1+p1*(x(1)+x(2)))); % clone 1

dxdt(2) = r2*x(2)*(1/(1+p2*(x(1)+x(2)))); % clone 2

dxdt(3) = r3*(x(1) + x (2)) - d3 * x(3); % M-protein

end