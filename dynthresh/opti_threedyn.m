%
% opti_threedyn.m
%
% function to optimize the RR of a three trial run with freedom to tune both
% thresholds.
%

function [h1m,h2m,h3m,h13s,RRdyn,RRconst] = opti_threedyn(eps,D,TD)

l = length(eps); h10=4; h20=3; h30=2; h130=2;
h1m = zeros(l,1); h2m = h1m; h3m=h1m; h13s = h1m; RRdyn=h1m; RRconst=h1m;
opts = optimset('Display','off','tolfun',1e-12);

for j=1:l, ep = eps(j);
    
    [X,fval] = fminsearch(@(X) RRfun(X,D,TD,ep),[h10; h20; h30],opts);
    h1m(j) = X(1); y20 = D*log(((1-ep)*exp(X(1)/D)+ep)/(ep*exp(X(1)/D)+(1-ep)));
    h2m(j) = max(X(2),y20); c2=1/(1+exp(-h2m(j)/D));
    y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
    h3m(j) = max(X(3),y30);
    h10=X(1)+0.01*rand; h20=X(2)+0.01*rand; h30=X(3)+0.01*rand;
    RRdyn(j) = -fval;
        
    [X,fval] = fminsearch(@(X) RRefun(X,D,TD,ep),h130,opts);
    h13s(j) = X; h130 = X; RRconst(j) = -fval;
    
end

function RR = RRfun(X,D,TD,ep)
h1=X(1); h2=X(2); h3=X(3);
c1 = 1/(1+exp(-h1/D)); DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));

if h2>y20
    c2 = 1/(1+exp(-h2/D));
    DT2 = h2*(1-exp(-h2/D))/(1+exp(-h2/D))-(1-2*ep)*(2*c1-1)*y20;
else
    c2 = (1-ep)*c1+ep*(1-c1); DT2 = 0;
end

y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));

if h3>y30
    c3 = 1/(1+exp(-h3/D));
    DT3 = h3*(1-exp(-h3/D))/(1+exp(-h3/D))-(1-2*ep)*(2*c2-1)*y30;
else
    c3 = (1-ep)*c2+ep*(1-c2); DT3 = 0;
end

RR = -(c1+c2+c3)/(DT1+DT2+DT3+3*TD);

function RR = RRefun(X,D,TD,ep)
h=X;
y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
RR = -3/((1-exp(-h/D))*(3*h-2*(1-2*ep)*y0)+3*TD*(1+exp(-h/D)));

