%
% sim_constthresh6.m
%
% simulate the two trial environment with constant threshold in trials 1-6
%

D = 1;
ep = 0;
TD = 2;
% plot_RRn;

hs = linspace(0,3,20);
hs = hs(end:-1:1);

Nsim = 1e5;
dt = 0.0005;
nos = sqrt(2*dt*D);
RRs = zeros(20,1);


for i=1:20, h = hs(i); cp1 = 0; cp2=0; cp3=0; cp4 = 0; cp5=0; cp6=0;
    DT1=0; DT2=0; DT3=0; DT4=0; DT5=0; DT6=0;
    y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
for j=1:Nsim
    
    % trial 1
    x=0; t=0;
    while abs(x)<h, x=x+dt+nos*randn; t=t+dt; end
    cp1=cp1+heaviside(x)/Nsim; DT1=DT1+t/Nsim;
    
    %trial 2
    x=sign(x)*y0; t=0; A=2*ceil(rand-ep)-1;
    while abs(x)<h, x=x+A*dt+nos*randn; t=t+dt; end
    cp2=cp2+heaviside(x*A)/Nsim; DT2=DT2+t/Nsim;
    
    %trial 3
    x=sign(x)*y0; t=0; A=(2*ceil(rand-ep)-1)*A;
    while abs(x)<h, x=x+A*dt+nos*randn; t=t+dt; end
    cp3=cp3+heaviside(x*A)/Nsim; DT3=DT3+t/Nsim;
    
    % trial 4
    x=sign(x)*y0; t=0; A=(2*ceil(rand-ep)-1)*A;
    while abs(x)<h, x=x+A*dt+nos*randn; t=t+dt; end
    cp4=cp4+heaviside(x*A)/Nsim; DT4=DT4+t/Nsim;
    
    %trial 5
    x=sign(x)*y0; t=0; A=(2*ceil(rand-ep)-1)*A;
    while abs(x)<h, x=x+A*dt+nos*randn; t=t+dt; end
    cp5=cp5+heaviside(x*A)/Nsim; DT5=DT5+t/Nsim;
    
    %trial 6
    x=sign(x)*y0; t=0; A=(2*ceil(rand-ep)-1)*A;
    while abs(x)<h, x=x+A*dt+nos*randn; t=t+dt; end
    cp6=cp6+heaviside(x*A)/Nsim; DT6=DT6+t/Nsim;
    
end
    RR = (cp1+cp2+cp3+cp4+cp5+cp6)/(DT1+DT2+DT3+DT4+DT5+DT6+6*TD);
    RRs(i) = RR;
    figure(1), hold on, plot(h,RR,'r.','markersize',40);
    pause(1e-9)
end