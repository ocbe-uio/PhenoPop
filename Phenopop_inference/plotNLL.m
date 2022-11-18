function plotNLL(DATA)
%plots negative log-likelihood and BIC versus population number.
Conc=[0 0.1 1 10];
Time=[0 24 48 72];

fval=zeros(6,1);
BIC=zeros(6,1);
LML=0.05;
[a,b,c]=size(DATA);
N=a*b*c;
for j=1:6
    j
    fval(j)=vectorized_inference_k_subpopulations_two_noise_levels(DATA,j,LML, Conc, Time);
    nparams=5*j+1;
    BIC(j)=nparams*log(N)+2*fval(j);
end

plot(1:6,fval,'-*','linewidth',5,'markersize',14)
xlabel('Number of Populations')
ylabel('Negative Log-likelihood')
ax = gca;
ax.FontSize = 15;
figure
plot(1:6,BIC,'-*','linewidth',5,'markersize',14)
xlabel('Number of Populations')
ylabel('BIC Score')    
ax = gca;
ax.FontSize = 15;