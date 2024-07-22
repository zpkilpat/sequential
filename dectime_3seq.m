%
% dectime_3seq.m
%
% plots the decision time as a function of eps, when there has been:
% RR, RA, AR, AA.

eps = linspace(0,0.5,200); el=length(eps);
TD = 2; % timedelay
D = 1;  % diff coef
TRR = zeros(el,1); TRA = zeros(el,1); TAR = zeros(el,1); TAA = zeros(el,1);

for j=1:el, ep = eps(j);
    RRifun = @(x)(1-exp(-x/D))*(x-(1-2*ep)*D*log(((1-ep)*exp(x/D)+ep)/...
        (ep*exp(x/D)+(1-ep))))+TD*(1+exp(-x/D));
    h = fminbnd(RRifun,0,10);
    y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
    Tp = h*(exp(h/D)+exp(-h/D)-2*exp(-y0/D))/(exp(h/D)-exp(-h/D))-y0;
    Tm = h*(exp(h/D)+exp(-h/D)-2*exp(+y0/D))/(exp(h/D)-exp(-h/D))+y0;
    cR = ((1-ep)*ep*exp(2*h/D)+(1-ep)*exp(h/D)+ep^2)/(1+exp(-h/D))/((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep));
    cA = ((1-ep)*ep*exp(2*h/D)+ep*exp(h/D)+(1-ep)^2)/(1+exp(-h/D))/((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep));
    TRR(j) = Tp*cR+Tm*(1-cR);
    TRA(j) = Tm*cR+Tp*(1-cR);
    TAR(j) = Tp*cA+Tm*(1-cA);
    TAA(j) = Tm*cA+Tp*(1-cA);
end

figure(1), hold on, plot(eps,TRR,'b','linewidth',8);
hold on, plot(eps,TRA,'r','linewidth',8);
hold on, plot(eps,TAR,'g','linewidth',8);
hold on, plot(eps,TAA,'m','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
% ylabel('$c_{RR}, c_{RA}, c_{AR}, c_{AA}$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); axis([0 0.5 0 1])