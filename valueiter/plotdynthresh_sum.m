%
% plotdynthresh_sum.m
%
% plots the optimal dynamic thresholds for a fixed ep as a function of the
% trial number.

D = 1;  % noise amplitude
TD = 1; % delay time
eps = 0.1;  % change rate
N = 15; % number of trials to run
cost = 0.25; % reward cost of time waiting
gam = 0.95;  % future rewards discounting

[Hm,Hc,rdyn,rconst] = opti_dynseq(eps,D,TD,cost,gam,N);

plot([1:N],Hc*ones(N,1),'linewidth',8);
hold on, plot([1:N],Hm,'linewidth',8);
set(gca,'fontsize',30);
xlabel('trial number','fontsize',33,'interpreter','latex')
ylabel('$\theta_{1:n}^{\rm max}$','fontsize',33,'interpreter','latex')