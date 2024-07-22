%
% constthreshmap_RR.m
%
% plots the reward rate as a function of a constant threshold.

D = 1;
TD = 2;
hs = linspace(0,2,1000);
eps = linspace(0,0.5,100);

[H,E] = meshgrid(hs,eps);

cp = 1./(1+exp(-H/D));
DT1 = H.*(1-exp(-H/D))./(1+exp(-H/D));
y0 = D*log(((1-E).*exp(H/D)+E)./(E.*exp(H/D)+(1-E)));
DT2 = (1-exp(-H/D)).*(H-(1-2*E).*y0)./(1+exp(-H/D));
RR = 2./((1-exp(-H/D)).*(2*H-(1-2*E).*y0)+2*TD*(1+exp(-H./D)));

hmax = zeros(100,1); RRmax = hmax;
for j=1:length(eps),
    [RRm,mi]=max(RR(j,:)); hmax(j) = hs(mi); RRmax(j)=RRm;
end

figure(1), hold on, pcolor(H,E,RR); colormap(hot); shading flat
hold on, plot(hmax,eps,'c','linewidth',5); colorbar
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$\epsilon$','fontsize',30,'interpreter','latex');
title('${\rm RR}_{1,2}(\theta, \epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);

figure(2), hold on, plot(eps,RRmax,'b','linewidth',8);
xlabel('$\epsilon$','fontsize',30,'interpreter','latex');
ylabel('${\rm RR}_{1,2}^{\rm max}(\epsilon)$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);
