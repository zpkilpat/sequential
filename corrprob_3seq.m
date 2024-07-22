%
% corrprob_3seq.m
%
% plots the correct probability as a function of eps, when there has been:
% RR, RA, AR, AA.

eps = linspace(0,0.5,200); el=length(eps);
TD = 2; % timedelay
D = 1;  % diff coef
cRR = zeros(el,1); cRA = zeros(el,1); cAR = zeros(el,1); cAA = zeros(el,1);

for j=1:el, ep = eps(j);
    RRifun = @(x)(1-exp(-x/D))*(x-(1-2*ep)*D*log(((1-ep)*exp(x/D)+ep)/...
        (ep*exp(x/D)+(1-ep))))+TD*(1+exp(-x/D));
    h = fminbnd(RRifun,0,10);
    y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
    pip = (1-ep)*exp(h/D)/((1-ep)*exp(h/D)+ep);
    pim = ep*exp(h/D)/(ep*exp(h/D)+1-ep);
    cR = ((1-ep)*ep*exp(2*h/D)+(1-ep)*exp(h/D)+ep^2)/(1+exp(-h/D))/((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep));
    cA = ((1-ep)*ep*exp(2*h/D)+ep*exp(h/D)+(1-ep)^2)/(1+exp(-h/D))/((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep));
    cRR(j) = pip*cR+pim*(1-cR);
    cRA(j) = pim*cR+pip*(1-cR);
    cAR(j) = pip*cA+pim*(1-cA);
    cAA(j) = pim*cA+pip*(1-cA);
end

figure(1), hold on, plot(eps,cRR,'b','linewidth',8);
hold on, plot(eps,cRA,'r','linewidth',8);
hold on, plot(eps,cAR,'g','linewidth',8);
hold on, plot(eps,cAA,'m','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('$c_{RR}, c_{RA}, c_{AR}, c_{AA}$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); axis([0 0.5 0 1])