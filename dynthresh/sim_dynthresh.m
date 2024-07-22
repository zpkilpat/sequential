%
% sim_dynthresh.m
%
% simulate the two trial environment with changing threshold in trial 1 and
% trial 2.
%

D = 1;
ep = 0.1;
TD = 2;
h1 = 0.5;
plot_dynthresh

h2s = linspace(0,2,20);
h2s = h2s(end:-1:1);

Nsim = 1e5;
dt = 0.0005;
nos = sqrt(2*dt*D);
cp1s = zeros(20,1); cp2s=cp1s; DT1s=cp1s; DT2s=cp1s;


for i=1:20, h2 = h2s(i); cp1 = 0; cp2=0; DT1=0; DT2=0;
    y0 = D*log(((1-ep)*exp(h1/D)+ep)./(ep*exp(h1/D)+(1-ep)));
for j=1:Nsim,
    
    % trial 1
    x=0; t=0;
    while abs(x)<h1, x=x+dt+nos*randn; t=t+dt; end
    cp1=cp1+heaviside(x)/Nsim; DT1=DT1+t/Nsim;
    
    %trial 2
    x=sign(x)*y0; t=0; A=2*ceil(rand-ep)-1;
    while abs(x)<h2, x=x+A*dt+nos*randn; t=t+dt; end
    cp2=cp2+heaviside(x*A)/Nsim; DT2=DT2+t/Nsim;
    
end
    RR = (cp1+cp2)/(DT1+DT2+2*TD);
    cp1s(i)=cp1; cp2s(i)=cp2; DT1s(i)=DT1; DT2s(i)=DT2; RRs(i) = RR;
%     figure(1), hold on, plot(h,DT1,'r.','markersize',40);
    figure(1), hold on, plot(h2,DT2,'m.','markersize',40);
%     figure(2), hold on, plot(h,cp1,'r.','markersize',40);
    figure(2), hold on, plot(h2,cp2,'m.','markersize',40);
    figure(3), hold on, plot(h2,(cp1+cp2)/(DT1+DT2+2*TD),'r.','markersize',40);
    pause(1e-9)
end