%
% opti_fouriterdyn.m
%
% function to optimize the RR of a four trial run with freedom to tune both
% thresholds. uses an iterative solver that sweeps multiple times along
% single parameter dimensions.
%

function [h1m,h2m,h3m,h4m,h14s,RRdyn,RRconst] = opti_fouriterdyn(eps,D,TD)

l = length(eps); h1=5; h2=4; h3=3; h4=3; h14=2;
h1m = zeros(l,1); h2m = h1m; h3m=h1m; h4m = h1m; h14s = h1m; RRdyn=h1m; RRconst=h1m;
opts = optimset('Display','off','tolfun',1e-15,'tolx',1e-15);

for j=1:l, ep = eps(j);
    
    for k=1:100
        h1 = fminbnd(@(X) RR1fun(X,h2,h3,h4,D,TD,ep),0,10,opts);
        y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
        h2 = fminbnd(@(X) RR2fun(X,h1,h3,h4,D,TD,ep),y20,10,opts);
        h2 = max(h2,y20); c2=1/(1+exp(-h2/D));
        y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
        h3 = fminbnd(@(X) RR3fun(X,h1,h2,h4,D,TD,ep),y30,10,opts);
        h3 = max(h3,y30); c3=1/(1+exp(-h3/D));
        y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
        [h4,fval] = fminbnd(@(X) RR4fun(X,h1,h2,h3,D,TD,ep),y40,10,opts);
    end
    h1m(j) = h1; y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
    h2m(j) = max(h2,y20); c2=1/(1+exp(-h2m(j)/D));
    y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
    h3m(j) = max(h3,y30); c3=1/(1+exp(-h3m(j)/D));
    y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
    h4m(j) = max(h4,y40);
    RRdyn(j) = -fval;
        
    [X,fval] = fminsearch(@(X) RRefun(X,D,TD,ep),h14,opts);
    h14s(j) = X; h14 = X; RRconst(j) = 1/fval;
    
end

function RR = RR1fun(X,h2,h3,h4,D,TD,ep)
h1=X;
c1 = 1/(1+exp(-h1/D)); DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
c2 = 1/(1+exp(-h2/D));
y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
DT2 = h2*(1-exp(-h2/D))/(1+exp(-h2/D))-(1-2*ep)*(2*c1-1)*y20;
c3 = 1/(1+exp(-h3/D));
y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
DT3 = h3*(1-exp(-h3/D))/(1+exp(-h3/D))-(1-2*ep)*(2*c2-1)*y30;
c4 = 1/(1+exp(-h4/D));
y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
DT4 = h4*(1-exp(-h4/D))/(1+exp(-h4/D))-(1-2*ep)*(2*c3-1)*y40;
% RR = (DT1+DT2+DT3+DT4+4*TD)/(c1+c2+c3+c4);
RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RR2fun(X,h1,h3,h4,D,TD,ep)
h2=X;
c1 = 1/(1+exp(-h1/D)); DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
c2 = 1/(1+exp(-h2/D));
y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
DT2 = h2*(1-exp(-h2/D))/(1+exp(-h2/D))-(1-2*ep)*(2*c1-1)*y20;
c3 = 1/(1+exp(-h3/D));
y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
DT3 = h3*(1-exp(-h3/D))/(1+exp(-h3/D))-(1-2*ep)*(2*c2-1)*y30;
c4 = 1/(1+exp(-h4/D));
y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
DT4 = h4*(1-exp(-h4/D))/(1+exp(-h4/D))-(1-2*ep)*(2*c3-1)*y40;
% RR = (DT1+DT2+DT3+DT4+4*TD)/(c1+c2+c3+c4);
RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RR3fun(X,h1,h2,h4,D,TD,ep)
h3=X;
c1 = 1/(1+exp(-h1/D)); DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
c2 = 1/(1+exp(-h2/D));
y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
DT2 = h2*(1-exp(-h2/D))/(1+exp(-h2/D))-(1-2*ep)*(2*c1-1)*y20;
c3 = 1/(1+exp(-h3/D));
y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
DT3 = h3*(1-exp(-h3/D))/(1+exp(-h3/D))-(1-2*ep)*(2*c2-1)*y30;
c4 = 1/(1+exp(-h4/D));
y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
DT4 = h4*(1-exp(-h4/D))/(1+exp(-h4/D))-(1-2*ep)*(2*c3-1)*y40;
% RR = (DT1+DT2+DT3+DT4+4*TD)/(c1+c2+c3+c4);
RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RR4fun(X,h1,h2,h3,D,TD,ep)
h4=X;
c1 = 1/(1+exp(-h1/D)); DT1 = h1*(1-exp(-h1/D))/(1+exp(-h1/D));
c2 = 1/(1+exp(-h2/D));
y20 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));
DT2 = h2*(1-exp(-h2/D))/(1+exp(-h2/D))-(1-2*ep)*(2*c1-1)*y20;
c3 = 1/(1+exp(-h3/D));
y30 = D*log(((1-ep)*c2+ep*(1-c2))/(ep*c2+(1-ep)*(1-c2)));
DT3 = h3*(1-exp(-h3/D))/(1+exp(-h3/D))-(1-2*ep)*(2*c2-1)*y30;
c4 = 1/(1+exp(-h4/D));
y40 = D*log(((1-ep)*c3+ep*(1-c3))/(ep*c3+(1-ep)*(1-c3)));
DT4 = h4*(1-exp(-h4/D))/(1+exp(-h4/D))-(1-2*ep)*(2*c3-1)*y40;
% RR = (DT1+DT2+DT3+DT4+4*TD)/(c1+c2+c3+c4);
RR = -(c1+c2+c3+c4)/(DT1+DT2+DT3+DT4+4*TD);

function RR = RRefun(X,D,TD,ep)
h=X;
y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
% RR = -4/((1-exp(-h/D))*(4*h-3*(1-2*ep)*y0)+4*TD*(1+exp(-h/D)));
RR = ((1-exp(-h/D))*(4*h-3*(1-2*ep)*y0)+4*TD*(1+exp(-h/D)))/4;

