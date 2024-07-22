%
% plotdynvtrial.m
%
% plots the optimal dynamic thresholds for a fixed ep as a function of the
% trial number.

D = 1;  % noise amplitude
TD = 2; % delay time
eps = linspace(0,0.5,100);
N = 10;
Hmtot = [];

for j=1:10
    [Hm,Hc,RRdyn,RRconst] = opti_Ndyn(eps,D,TD,j);
    Hmtot = [Hmtot;Hm];
end

ei = 50; cc = Hmtot(:,ei); ind=1;
for j=1:10, ind=j*(j-1)/2+1;
    hold on, plot([1:j],cc(ind:ind+j-1),'color',[1/2 j/10 1/2],'linewidth',8);
end
set(gca,'fontsize',30,'xtick',[1:10]);
xlabel('trial number','fontsize',33,'interpreter','latex')
ylabel('$\theta_{1:n}^{\rm max}$','fontsize',33,'interpreter','latex')