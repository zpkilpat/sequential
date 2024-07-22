%
% opti_twodyn.m
%
% function to optimize the RR of a two trial run with freedom to tune both
% thresholds.
%

function [h1m,h2m,h12s,RRdyn] = opti_twodyn(eps,D,TD)

l = length(eps); h10=1; h20=1; h120=1;
h1m = zeros(l,1); h2m = h1m; h12s = h1m; RRdyn=h1m;
opts = optimset('Display','off','tolfun',1e-9,'algorithm','levenberg-marquardt');

for j=1:l, ep = eps(j);
    
    [X,fval] = fminsearch(@(X) RRfun(X,D,TD,ep),[h10; h20],opts);
    h1m(j) = X(1); y20 = D*log(((1-ep)*exp(X(1)/D)+ep)/(ep*exp(X(1)/D)+(1-ep)));
    h2m(j) = max(X(2),y20); h10=X(1); h20=X(2);
    RRdyn(j) = -fval;
        
    X = fminsearch(@(X) RRefun(X,D,TD,ep),h120,opts);
    h12s(j) = X; h120 = X;
    
end

function RR = RRfun(X,D,TD,ep)
h1=X(1); h2=X(2);
Ep1=1+exp(-h1/D); Ep2=1+exp(-h2/D); Em1=1-exp(-h1/D); Em2=1-exp(-h2/D);
y0 = D*log(((1-ep)*exp(h1/D)+ep)/(ep*exp(h1/D)+(1-ep)));

if h2>y0
    RR = -(Ep1+Ep2)/(h1*Ep2*Em1+h2*Ep1*Em2+2*TD*Ep1*Ep2-(1-2*ep)*Ep2*Em1*y0);
else
    RR = -(2-ep*(1-exp(-h1/D)))/(h1*(1-exp(-h1/D))+2*TD*(1+exp(-h1/D)));
end

function RR = RRefun(X,D,TD,ep)
h=X;
y0 = D*log(((1-ep)*exp(h/D)+ep)/(ep*exp(h/D)+(1-ep)));
RR = -2/((1-exp(-h/D))*(2*h-(1-2*ep)*y0)+2*TD*(1+exp(-h/D)));

