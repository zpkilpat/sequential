%
% opti_fivedyn.m
%
% function to optimize the RR of a five trial run with freedom to tune both
% thresholds.
%

function [h1m,h2m,h3m,h4m,h5m,h15s,RRdyn,RRconst] = opti_fivedyn(eps,D,TD)

l = length(eps); h10=5; h20=4; h30=3; h40=3; h50=3; h150=2;
h1m = zeros(l,1); h2m = h1m; h3m=h1m; h4m = h1m; h5m = h1m; h15s = h1m; RRdyn=h1m; RRconst=h1m;
opts = optimset('Display','off','tolfun',1e-15,'tolx',1e-15);

for j=1:l, ep = eps(j);
    
    [X,fval] = fminsearch(@(X) RRfun(X,D,TD,ep),[h10; h20; h30; h40; h50],opts);
    h1m(j) = X(1); y20 = D*log(((1-ep)*exp(X(1)/D)+ep)/(ep*exp(X(1)/D)+(1-ep)));
    h2m(j) = max(X(2),y20); c2=1/(1+exp(-h2m(j)/D));
    y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
    h3m(j) = max(X(3),y30); c3=1/(1+exp(-h3m(j)/D));
    y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
    h4m(j) = max(X(4),y40); c4=1/(1+exp(-h4m(j)/D));
    y50 = D*log(((1-ep)*c4+ep*(1-c4))/(ep*c4+(1-ep)*(1-c4)));
    h5m(j) = max(X(5),y50);
    h10=h1m(j); h20=h2m(j); h30=h3m(j); h40=h4m(j); h50=h5m(j);
    RRdyn(j) = 1/fval;
        
    [X,fval] = fminsearch(@(X) RRefun(X,D,TD,ep),h150,opts);
    h15s(j) = X; h150 = X; RRconst(j) = 1/fval;
    
end

function RR = RRfun(X,D,TD,ep)
h1=X(1); h2=X(2); h3=X(3); h4=X(4); h5=X(5);
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

y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));

if h4>y40
    c4 = 1/(1+exp(-h4/D));
    DT4 = h4*(1-exp(-h4/D))/(1+exp(-h4/D))-(1-2*ep)*(2*c3-1)*y40;
else
    c4 = (1-ep)*c3+ep*(1-c3); DT4 = 0;
end

y50 = D*log(((1-ep)*c4+ep*(1-c4))/(ep*c4+(1-ep)*(1-c4)));

if h5>y50
    c5 = 1/(1+exp(-h5/D));
    DT5 = h5*(1-exp(-h5/D))/(1+exp(-h5/D))-(1-2*ep)*(2*c4-1)*y50;
else
    c5 = (1-ep)*c4+ep*(1-c4); DT5 = 0;
end

RR = (DT1+DT2+DT3+DT4+DT5+5*TD)/(c1+c2+c3+c4+c5);
% RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RRefun(X,D,TD,ep)
h=X;
y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
% RR = -4/((1-exp(-h/D))*(4*h-3*(1-2*ep)*y0)+4*TD*(1+exp(-h/D)));
RR = ((1-exp(-h/D))*(5*h-4*(1-2*ep)*y0)+5*TD*(1+exp(-h/D)))/5;

