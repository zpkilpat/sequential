%
% max_Ndyn.m
%
% script for maximizing the RR for an N trial dynamic threshold system.
%

D = 1;  % noise amplitude
TD = 2; % delay time
eps = linspace(0,0.5,100); eps = eps(end:-1:1);
N = 10; % number of trials in a row in an experiment

% run the function to return the max's for a range of eps
[Hm,Hc,RRdyn,RRconst] = opti_Ndyn(eps,D,TD,N);

figure(1), hold on, plot(eps,Hc,'k','linewidth',3);
for j=1:N, hold on, plot(eps,Hm(j,:),'linewidth',8); end
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('$\theta_{1:N}^{\rm max}$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);
set(gca,'xtick',0:0.1:0.5)

figure(2), hold on, plot(eps,RRconst(1)+0*eps,'k','linewidth',3);
hold on, plot(eps,RRconst,'k','linewidth',3);
hold on, plot(eps,RRdyn,'b','linewidth',8);
set(gca,'xtick',0:0.1:0.5)
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('${\rm RR}_{1:N}^{\rm max}(\epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);