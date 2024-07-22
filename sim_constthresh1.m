%
% sim_constthresh1.m
%
% simulate the two trial environment with constant threshold in trial 1.
%

D = 1;
TD = 2;
% constthresh_RR;

hs = linspace(0,3,20);
hs = hs(end:-1:1);

Nsim = 1e5;
dt = 0.0005;
nos = sqrt(2*dt*D);
RRs = zeros(20,1);


for i=1:20, h = hs(i); cp1 = 0; DT1=0;
for j=1:Nsim
    
    % trial 1
    x=0; t=0;
    while abs(x)<h, x=x+dt+nos*randn; t=t+dt; end
    cp1=cp1+heaviside(x)/Nsim; DT1=DT1+t/Nsim;
    
end
    RR = cp1/(DT1+TD); RRs(i)=RR;
    figure(1), hold on, plot(h,RR,'r.','markersize',40);
    pause(1e-9)
end