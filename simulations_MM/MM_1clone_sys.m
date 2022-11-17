function dxdt = MM_1clone_sys(t,x,teta)
%ode system for MM dynamics with 2 clones

%notation of model parameters   
r1=teta(1);
p1=teta(2);
r3=teta(3);
d3=teta(4);


dxdt=zeros(2,1);

dxdt(1) = r1*x(1)*(1/(1+p1*(x(1)))); % clone 1

dxdt(2) = r3*x(1) - d3 * x(2); % M-protein

end