%
% dynthresh_RR.m
%
% plots the matrix of reward rates in the case of a dynamic threshold.
%

D = 1;
TD = 2;
ep = 0.3;

h1 = linspace(0,2,1000);
h2 = linspace(0,2,1000);
[H1,H2] = meshgrid(h1,h2);

Ep1=1+exp(-H1/D); Ep2=1+exp(-H2/D); Em1=1-exp(-H1/D); Em2=1-exp(-H2/D);

y0 = D*log(((1-ep)*exp(H1/D)+ep)./(ep*exp(H1/D)+(1-ep)));
y0vec = D*log(((1-ep)*exp(h1/D)+ep)./(ep*exp(h1/D)+(1-ep)));

RR = (Ep1+Ep2)./(H1.*Ep2.*Em1+H2.*Ep1.*Em2+2*TD*Ep1.*Ep2-(1-2*ep)*Ep2.*Em1.*y0);
RRs = (2-ep*(1-exp(-H1/D)))./(H1.*(1-exp(-H1/D))+2*TD*(1+exp(-H1/D)));
RR(H2<y0) = RRs(H2<y0);
RReq = diag(RR);

[RRmax,maxin] = max(RR(:));
h1max=H1(maxin); h2max=H2(maxin);
if h2max==0, h2max=D*log(((1-ep)*exp(h1max/D)+ep)./(ep*exp(h1max/D)+(1-ep))); end
[junk,maxin]=max(RReq); h12eq = h1(maxin);

pcolor(H1,H2,RR); shading flat, colormap(hot), colorbar
hold on, plot(h1,h1,'k','linewidth',2)
hold on, plot(h1,y0vec,'b--','linewidth',5);
hold on, plot(h1max,h2max,'r.','markersize',50);
%hold on, plot(h12eq,h12eq,'ro','linewidth',4,'markersize',15);
title('${\rm RR}_{1,2}(\theta_1,\theta_2;\epsilon=0)$','fontsize',30,'interpreter','latex');
xlabel('$\theta_1$','interpreter','latex','fontsize',30);
ylabel('$\theta_2$','interpreter','latex','fontsize',30);
set(gca,'fontsize',30);