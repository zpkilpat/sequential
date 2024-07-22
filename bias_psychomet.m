%
% bias_psychomet.m
%
% plots the psychometric function as a function of coherence (1/D),
% choosing for each coherence the optimal \theta for an infinite number of
% trials.

TD = 2; % time delay
ep = 0.25;   % switch rate
Dis = linspace(0.01,10,200);
Ds = 1./Dis;
cp0 = zeros(200,1); cpp=cp0; cpm=cp0;

for j=1:200, D=Ds(j);
    RRifun = @(x)(1-exp(-x/D))*(x-(1-2*ep)*D*log(((1-ep)*exp(x/D)+ep)/...
        (ep*exp(x/D)+(1-ep))))+TD*(1+exp(-x/D));
    h = fminbnd(RRifun,0,10);
    y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
    cp0(200+j) = 1/(1+exp(-h/D));
    cpp(200+j) = (1-exp(-(y0+h)/D))/(1-exp(-2*h/D));
    cpm(200+j) = (1-exp(-(-y0+h)/D))/(1-exp(-2*h/D));
end
cp0(1:200)=1-cp0(400:-1:201);
cpp(1:200)=1-cpm(400:-1:201);
cpm(1:200)=1-cpp(400:-1:201);

% figure(1), hold on, plot([0 0],[0 1],'k','linewidth',2);
% hold on, plot([-Dis(200:-1:1), Dis],cp0,'k','linewidth',8);
% hold on, plot([-Dis(200:-1:1), Dis],cpp,'b','linewidth',8);
% hold on, plot([-Dis(200:-1:1), Dis],cpm,'b','linewidth',8);
% xlabel('coherence ($g/D$)','fontsize',30,'interpreter','latex');
% ylabel('P($d_j = +1$) ','fontsize',30,'interpreter','latex')
% set(gca,'fontsize',30);