%
% plot_RRn.m
%
% plot the reward rate function for n trials for constant threshold.
%

% n = 6; % no. of trials
D = 1; % amp. of noise
ep = 0;  % change rate
TD = 2;     % delay time
hs = linspace(0,3,1000);

y0 = D*log(((1-ep)*exp(hs/D)+ep)./(ep*exp(hs/D)+(1-ep)));
RR = n./((1-exp(-hs/D)).*(n*hs-(n-1)*(1-2*ep)*y0)+n*TD*(1+exp(-hs/D)));
RRinf = 1./((1-exp(-hs/D)).*(hs-(1-2*ep)*y0)+TD*(1+exp(-hs/D)));

figure(1), hold on, plot(hs,RR,'b','linewidth',5);
% hold on, plot(hs,RRinf,'b','linewidth',5);
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('${\rm RR}_{1:n}(\theta; \epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);