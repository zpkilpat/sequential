%
% epsassume_colorplot.m
%
% make a colormap plot that shows the RR as a function of eps_assume and
% eps_true.
%

D = 1;
TD = 2;
eas = linspace(0,0.5,100);
ets = linspace(0,0.5,100);
RRm = zeros(100);
hs = linspace(0,5,1e5);
mi = zeros(100,1);

for j=1:100, et = ets(j); j
    for k=1:100, ea = eas(k);
c1 = 1./(1+exp(-hs/D));
DT1 = hs.*(1-exp(-hs/D))./(1+exp(-hs/D));

pp = (1-et)*c1+et*(1-c1); pm = et*c1+(1-et)*(1-c1);
y0 = D*log(((1-ea)*exp(hs/D)+ea)./(ea*exp(hs/D)+(1-ea)));

cp = (1-exp(-(y0+hs)/D))./(1-exp(-2*hs/D));
cm = (1-exp(-(-y0+hs)/D))./(1-exp(-2*hs/D));
Tp = hs.*(exp(hs/D)+exp(-hs/D)-2*exp(-y0/D))./(exp(hs/D)-exp(-hs/D))-y0;
Tm = hs.*(exp(hs/D)+exp(-hs/D)-2*exp(+y0/D))./(exp(hs/D)-exp(-hs/D))+y0;

c2 = pp.*cp+pm.*cm; DT2 = pp.*Tp+pm.*Tm;
RR = (c1+c2)./(DT1+DT2+2*TD);

RRs(j,k) = max(RR);

    end
    plot(eas,RRs(j,:)), hold on, pause(0.001);
    [junk,mi(j)] = max(RRs(j,:));
end
close

[ET,EA] = meshgrid(ets,eas);

figure(1), pcolor(ET,EA,RRs), shading flat, colormap(hot);
hold on, plot(eas(mi),ets,'linewidth',4);

% figure(1),
% hold on, plot(hs,DT1,'b','linewidth',5);
% hold on, plot(hs,DT2,'c','linewidth',5);
% xlabel('$\theta$','fontsize',30,'interpreter','latex');
% ylabel('$DT_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
% set(gca,'fontsize',30);
% 
% figure(2), hold on, plot(hs,c1,'b','linewidth',5);
% hold on, plot(hs,c2,'c','linewidth',5);
% xlabel('$\theta$','fontsize',30,'interpreter','latex');
% ylabel('$c_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
% set(gca,'fontsize',30);
% 
% figure(3), hold on, plot(hs,RR,'b','linewidth',5);
% xlabel('$\theta$','fontsize',30,'interpreter','latex');
% ylabel('$RR_{1:2}(\theta)$','fontsize',30,'interpreter','latex');
% set(gca,'fontsize',30);