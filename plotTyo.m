%
% plotTyo.m
%
% shows the dependence of Tyo on epsilon.

D = 0.1;
h = 0.1;
eps = linspace(0,0.5);

y0 = D*log(((1-eps)*exp(h/D)+eps)./(eps*exp(h/D)+(1-eps)));

Ty0p = h*coth(h/D)-h*csch(h/D)*exp(-y0/D)-y0;
Ty0m = h*coth(h/D)-h*csch(h/D)*exp(+y0/D)+y0;
dif = 2*h*csch(h/D)*sinh(y0/D)-2*y0;
dif2 = 2*h*csch(h/D)*sinh(y0/D);

figure(1), hold on, plot(eps,Ty0p,'k','linewidth',4);
hold on, plot(eps,Ty0m,'b','linewidth',4);

figure(2), hold on, plot(eps,Ty0p-Ty0m,'k','linewidth',4);
hold on, plot(eps,dif,'r--','linewidth',4);

figure(3), hold on, plot(eps,dif2,'g','linewidth',4);

x = linspace(0,1);
f = x-tanh(x);
figure(4), hold on, plot(x,f,'b','linewidth',4);
