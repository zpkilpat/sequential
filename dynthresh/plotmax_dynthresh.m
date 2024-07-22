%
% plotmax_dynthresh.m
%
% plots the matrix of reward rates in the case of a dynamic threshold.
%

D = 1;
TD = 2;
eps = linspace(0,0.5,100);
h1m = zeros(100,1); h2m = h1m; h12s = h1m; RRdyn=h1m; RRconst=h1m;

for j=1:100, ep = eps(j);
h1 = linspace(0,2,2000);
h2 = linspace(0,2,2000);
[H1,H2] = meshgrid(h1,h2);

Ep1=1+exp(-H1/D); Ep2=1+exp(-H2/D); Em1=1-exp(-H1/D); Em2=1-exp(-H2/D);

y0 = D*log(((1-ep)*exp(H1/D)+ep)./(ep*exp(H1/D)+(1-ep)));
y0vec = D*log(((1-ep)*exp(h1/D)+ep)./(ep*exp(h1/D)+(1-ep)));

RR = (Ep1+Ep2)./(H1.*Ep2.*Em1+H2.*Ep1.*Em2+2*TD*Ep1.*Ep2-(1-2*ep)*Ep2.*Em1.*y0);
RRs = (2-ep*(1-exp(-H1/D)))./(H1.*(1-exp(-H1/D))+2*TD*(1+exp(-H1/D)));
RR(H2<y0) = RRs(H2<y0);
RReq = diag(RR);

[RRmax,maxin] = max(RR(:)); RRdyn(j) = RRmax;
h1max=H1(maxin); h2max=H2(maxin);
if h2max==0, h2max=D*log(((1-ep)*exp(h1max/D)+ep)./(ep*exp(h1max/D)+(1-ep))); end
[RR12,maxin]=max(RReq); h12eq = h1(maxin); RRconst(j)=RR12;

h1m(j) = h1max; h2m(j) = h2max; h12s(j) = h12eq;

end

figure(1), hold on, plot(eps,h12s,'k','linewidth',3);
hold on, plot(eps,h1m,'b','linewidth',8);
hold on, plot(eps,h2m,'r','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('$\theta_{1,2}^{\rm max}$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

figure(2), hold on, plot(eps,RRdyn,'b','linewidth',8);
hold on, plot(eps,RRconst,'b','linewidth',3);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('${\rm RR}_{1,2}^{\rm max}(\epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

