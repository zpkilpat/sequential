%
% max_fivedyn.m
%
% script for maximizing the RR for a five trial dynamic threshold system.
%

D = 1;  % noise amplitude
TD = 2; % delay time
eps = linspace(0,0.5,100);
% eps = eps(end:-1:1);

% run the function to return the max's for a range of eps
[h1m,h2m,h3m,h4m,h5m,h15s,RRdyn,RRconst] = opti_fivedyn(eps,D,TD);

figure(1), hold on, plot(eps,h15s,'k','linewidth',3);
hold on, plot(eps,h1m,'b','linewidth',8);
hold on, plot(eps,h2m,'r','linewidth',8);
hold on, plot(eps,h3m,'g','linewidth',8);
hold on, plot(eps,h4m,'c','linewidth',8);
hold on, plot(eps,h5m,'m','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('$\theta_{1:5}^{\rm max}$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

figure(2), hold on, plot(eps,RRconst,'k','linewidth',3);
hold on, plot(eps,RRdyn,'b','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('${\rm RR}_{1:5}^{\rm max}(\epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);