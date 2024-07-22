%
% sim_dyn3RRmax.m
%
% run simulations to generate estimates of RR at predicted maximum values


Nsim = 1e5;
dt = 0.0005;
nos = sqrt(2*dt*D);
cp1s = zeros(20,1); cp2s=cp1s; DT1s=cp1s; DT2s=cp1s;


for i=1:20, ep = epsub(i); h1=h1sub(i); h2=h2sub(i); h3=h3sub(i);
    cp1 = 0; cp2=0; cp3=0; DT1=0; DT2=0; DT3=0;
    y20 = D*log(((1-ep)*exp(h1/D)+ep)./(ep*exp(h1/D)+(1-ep)));
    c1 = 1/(1+exp(-h1/D));
    if h2>y20, c2 = 1/(1+exp(-h2/D)); else, c2 = (1-ep)*c1+ep*(1-c1); end
    y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
    
for j=1:Nsim,
    
    % trial 1
    x=0; t=0;
    while abs(x)<h1, x=x+dt+nos*randn; t=t+dt; end
    cp1=cp1+heaviside(x)/Nsim; DT1=DT1+t/Nsim;
    
    %trial 2
    x=sign(x)*y20; t=0; A=2*ceil(rand-ep)-1;
    while abs(x)<h2, x=x+A*dt+nos*randn; t=t+dt; end
    cp2=cp2+heaviside(x*A)/Nsim; DT2=DT2+t/Nsim;
    
    %trial 3
    x=sign(x)*y30; t=0; A=2*ceil(rand-ep)-1;
    while abs(x)<h3, x=x+A*dt+nos*randn; t=t+dt; end
    cp3=cp3+heaviside(x*A)/Nsim; DT3=DT3+t/Nsim;
    
end
    RR = (cp1+cp2+cp3)/(DT1+DT2+DT3+3*TD);
    cp1s(i)=cp1; cp2s(i)=cp2; cp3s(i)=cp3; DT1s(i)=DT1; DT2s(i)=DT2; DT3s(i)=DT3; RRs(i) = RR;
    figure(2), hold on, plot(epsub(i),RR,'r.','markersize',40);
    pause(1e-9)
end