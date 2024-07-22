%
% opti_Ndyn.m
%
% function to optimize the RR of a five trial run with freedom to tune both
% thresholds.
%

function [Hm,Hc,RRdyn,RRconst] = opti_Ndyn(eps,D,TD,N)

l = length(eps); h0=0.795*ones(N,1); hc0=2;
Hm = zeros(N,l); Hc = zeros(1,l); RRdyn=Hc; RRconst=Hc;
opts = optimset('Display','off','tolfun',1e-15,'tolx',1e-15);

for j=1:l, ep = eps(j);
    
    [X,fval] = fminsearch(@(X) RRfun(X,D,TD,ep,N),h0,opts);
    Hm(1,j) = X(1); c=1/(1+exp(-Hm(1,j)/D));
    for k=2:N
    y0 = D*log(((1-ep)*c+ep*(1-c))/(ep*c+(1-ep)*(1-c)));
    Hm(k,j) = max(X(k),y0); c=1/(1+exp(-Hm(k,j)/D));
    end
    h0=Hm(:,j); RRdyn(j) = 1/fval;
        
    [X,fval] = fminsearch(@(X) RRefun(X,D,TD,ep,N),hc0,opts);
    Hc(j) = X; hc0 = X; RRconst(j) = 1/fval;
    
end

function RR = RRfun(X,D,TD,ep,N)
h=X; c=zeros(N,1); DT=zeros(N,1);
c(1) = 1/(1+exp(-h(1)/D)); DT(1) = h(1)*(1-exp(-h(1)/D))/(1+exp(-h(1)/D));

for k=2:N
y0 = D*log(((1-ep)*c(k-1)+ep*(1-c(k-1)))/(ep*c(k-1)+(1-ep)*(1-c(k-1))));
if h(k)>y0
    c(k) = 1/(1+exp(-h(k)/D));
    DT(k) = h(k)*(1-exp(-h(k)/D))/(1+exp(-h(k)/D))-(1-2*ep)*(2*c(k-1)-1)*y0;
else
    c(k) = (1-ep)*c(k-1)+ep*(1-c(k-1)); DT(k) = 0;
end
end
RR = (sum(DT)+N*TD)/sum(c);
% RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RRefun(X,D,TD,ep,N)
h=X;
y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
% RR = -4/((1-exp(-h/D))*(4*h-3*(1-2*ep)*y0)+4*TD*(1+exp(-h/D)));
RR = ((1-exp(-h/D))*(N*h-(N-1)*(1-2*ep)*y0)+N*TD*(1+exp(-h/D)))/N;

