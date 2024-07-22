%
% constthresh_RR.m
%
% plots the reward rate as a function of a constant threshold.

% D = 1;
% ep = 0.1;
% TD = 2;
hs = linspace(0,2,1000);

cp = 1./(1+exp(-hs/D));
DT1 = hs.*(1-exp(-hs/D))./(1+exp(-hs/D));
y0 = D*log(((1-ep)*exp(hs/D)+ep)./(ep*exp(hs/D)+(1-ep)));
DT2 = (1-exp(-hs/D)).*(hs-(1-2*ep)*y0)./(1+exp(-hs/D));
RR = 2./((1-exp(-hs/D)).*(2*hs-(1-2*ep)*y0)+2*TD*(1+exp(-hs./D)));
% RRo = 2*cp./(DT1+DT2+2*TD);

figure(1),
% hold on, plot(hs,DT1,'b','linewidth',5);
hold on, plot(hs,DT2,'c','linewidth',5);
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$DT_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

% figure(2), hold on, plot(hs,cp,'b','linewidth',5);
% xlabel('$\theta$','fontsize',30,'interpreter','latex');
% ylabel('$c_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
% set(gca,'fontsize',30);

figure(3), hold on, plot(hs,RR,'b','linewidth',5);
% hold on, plot(hs,RRo,'r','linewidth',5);
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$RR_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);
