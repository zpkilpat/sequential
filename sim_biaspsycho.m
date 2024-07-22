%
% sim_biaspsycho.m
%
% plots the psychometric function as a function of coherence (1/D),
% choosing for each coherence the optimal \theta for an infinite number of
% trials. plots are generated from simulation.

TD = 2; % time delay
ep = 0.25;   % switch rate
Dis = linspace(-9.5,9.5,20);
Ds = 1./abs(Dis);
gs = [-ones(1,10), ones(1,10)];
cp0 = zeros(20,1); cpp=cp0; cpm=cp0;
Nsim = 1e5; dt=0.0005;

for j=1:20, D=Ds(j); nos=sqrt(2*dt*D); g=gs(j);
    RRifun = @(x)(1-exp(-x/D))*(x-(1-2*ep)*D*log(((1-ep)*exp(x/D)+ep)/...
        (ep*exp(x/D)+(1-ep))))+TD*(1+exp(-x/D));
    h = fminbnd(RRifun,0,10);
    y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
    
    for k=1:Nsim
    % unbiased trials
    x=0; t=0;
    while abs(x)<h, x=x+g*dt+nos*randn; t=t+dt; end
    cp0(j)=cp0(j)+heaviside(x)/Nsim;
    
    % up-biased trials
    x=y0; t=0;
    while abs(x)<h, x=x+g*dt+nos*randn; t=t+dt; end
    cpp(j)=cpp(j)+heaviside(x)/Nsim;
    
    % down-biased trials
    x=-y0; t=0;
    while abs(x)<h, x=x+g*dt+nos*randn; t=t+dt; end
    cpm(j)=cpm(j)+heaviside(x)/Nsim;
    end
    
end

figure(1), hold on, plot(Dis,cp0,'r.','markersize',40);
hold on, plot(Dis,cpp,'r.','markersize',40);
hold on, plot(Dis,cpm,'r.','markersize',40);