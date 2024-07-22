%
% opti_dynseq.m
%
% function to optimize the RR of a five trial run with freedom to tune both
% thresholds.
%

function [Hm,Hc,rdyn,rconst] = opti_dynseq(eps,D,TD,cost,gam,N)

l = length(eps); h0=0.795*ones(N,1); hc0=2;
Hm = zeros(N,l); Hc = zeros(1,l); rdyn=Hc; rconst=Hc;
opts = optimset('Display','off','tolfun',1e-21,'tolx',1e-21);

for j=1:l, ep = eps(j);
    
    [X,fval] = fminsearch(@(X) rfun(X,D,TD,ep,cost,gam,N),h0,opts);
    Hm(1,j) = X(1); c=1/(1+exp(-Hm(1,j)/D));
    for k=2:N
    y0 = D*log(((1-ep)*c+ep*(1-c))/(ep*c+(1-ep)*(1-c)));
    Hm(k,j) = max(X(k),y0); c=1/(1+exp(-Hm(k,j)/D));
    end
    h0=Hm(:,j); rdyn(j) = -fval;
        
    [X,fval] = fminsearch(@(X) rcfun(X,D,TD,ep,cost,gam),hc0,opts);
    Hc(j) = X; hc0 = X; rconst(j) = -fval;
    
end

% compute negative of the objective function (so that we can minimize
% rather than maximize). this is the general case of dynamic thresholds.
function rval = rfun(X,D,TD,ep,cost,gam,N)
h=X; c=zeros(N,1); DT=zeros(N,1);
c(1) = 1/(1+exp(-h(1)/D)); DT(1) = h(1)*(1-exp(-h(1)/D))/(1+exp(-h(1)/D))+TD;
rval = -c(1)+cost*DT(1);

for k=2:N
y0 = D*log(((1-ep)*c(k-1)+ep*(1-c(k-1)))/(ep*c(k-1)+(1-ep)*(1-c(k-1))));
if h(k)>y0
    c(k) = 1/(1+exp(-h(k)/D));
    DT(k) = h(k)*(1-exp(-h(k)/D))/(1+exp(-h(k)/D))-(1-2*ep)*(2*c(k-1)-1)*y0+TD;
else
    c(k) = (1-ep)*c(k-1)+ep*(1-c(k-1)); DT(k) = TD;
end
rval=rval+gam^(k-1)*(cost*DT(k)-c(k));
end

% compute negative of the objective function when fixing threshold constant
% across trials.
function rval = rcfun(X,D,TD,ep,cost,gam)
h=X;
c = 1/(1+exp(-X/D));
y0 = D*log(((1-ep)*exp(X/D)+ep)/(ep*exp(X/D)+(1-ep)));
rval = cost*(X*(1-exp(-X/D))/(1+exp(-X/D))+TD)-c+gam/(1-gam)*...
    (cost*((X-(1-2*ep)*y0)*(1-exp(-X/D))*c+TD)-c);

