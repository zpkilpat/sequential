%
% plot_dynthresh.m
%
% plots the decision time and reward rate for fixed theta_1 as a function
% of theta_2

% D = 1;
% ep = 0.5;
% TD = 2;
% h1 = 0.5;
hs = linspace(0,2,1000);

c1 = 1./(1+exp(-h1/D));
DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
c2 = 1./(1+exp(-hs/D));
y0 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
c2(hs<y0) = (1-ep)*c1+ep*(1-c1);
DT2 = hs.*(1-exp(-hs/D))./(1+exp(-hs/D))-(1-2*ep)*(1-exp(-h1/D))*y0/(1+exp(-h1/D));
DT2(hs<y0) = 0;
RR = (c1+c2)./(DT1+DT2+2*TD);

figure(1), hold on, plot(hs,DT2,'b','linewidth',5);
xlabel('$\theta_2$','fontsize',30,'interpreter','latex');
ylabel('$DT_{2}(\theta)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

figure(2), hold on, plot(hs,c2,'b','linewidth',5);
xlabel('$\theta_2$','fontsize',30,'interpreter','latex');
ylabel('$c_{2}(\theta)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

figure(3), hold on, plot(hs,RR,'b','linewidth',5);
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$RR_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);