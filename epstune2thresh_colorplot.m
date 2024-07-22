%
% epstune2thresh_colorplot.m
%
% make a colormap plot that shows the RR as a function of eps_assume and
% eps_true.
%

D = 1;
TD = 2;
eas = linspace(0,0.5,100);
ets = linspace(0,0.5,100);
RRm = zeros(100); mi = zeros(100,1);
h1 = linspace(0,2,1000); h2 = linspace(0,2,1000);  [H1,H2] = meshgrid(h1,h2);


for j=1:100, et = ets(j); j
    for k=1:100, ea = eas(k);
Ep1=1+exp(-H1/D); Ep2=1+exp(-H2/D);
c1 = 1./Ep1; DT1 = H1.*Em1./Ep1;

pp = (1-et)*c1+et*(1-c1); pm = et*c1+(1-et)*(1-c1);
y0 = D*log(((1-ea)*exp(H1/D)+ea)./(ea*exp(H1/D)+(1-ea)));

cp = (1-exp(-(y0+H2)/D))./(1-exp(-2*H2/D));
cm = (1-exp(-(-y0+H2)/D))./(1-exp(-2*H2/D));
Tp = H2.*(exp(H2/D)+exp(-H2/D)-2*exp(-y0/D))./(exp(H2/D)-exp(-H2/D))-y0;
Tm = H2.*(exp(H2/D)+exp(-H2/D)-2*exp(+y0/D))./(exp(H2/D)-exp(-H2/D))+y0;

c2 = pp.*cp+pm.*cm; DT2 = pp.*Tp+pm.*Tm;
c2(H2<y0)=(1-et)*c1(H2<y0)+et*(1-c1(H2<y0)); DT2(H2<y0)=0;
RR = (c1+c2)./(DT1+DT2+2*TD);

RRs(j,k) = max(RR(:));

    end
    plot(eas,RRs(j,:)), hold on, pause(0.001);
    [junk,mi(j)] = max(RRs(j,:));
end
close

[ET,EA] = meshgrid(ets,eas);

figure(1), pcolor(ET,EA,RRs'), shading flat, colormap(hot);
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