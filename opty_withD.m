%
% opty_withD.m
%
% plots the optimal initial condition as a function of D

ep = 0.4;
TD = 2; % timedelay
Ds = linspace(0.0001,5);  % diff coef
y0s = zeros(100,1);

for j=2:100, D = Ds(j);
    RRifun = @(x)(1-exp(-x/D))*(x-(1-2*ep)*D*log(((1-ep)*exp(x/D)+ep)/...
        (ep*exp(x/D)+(1-ep))))+TD*(1+exp(-x/D));
    h = fminbnd(RRifun,0,10);
    y0s(j) = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
end

figure(1), hold on, plot(Ds,y0s,'linewidth',8);
xlabel('$D$','fontsize',30,'interpreter','latex');
ylabel('$y_0$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);