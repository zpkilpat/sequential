%
% simDDM.m
%
% simulate a single DDM with a tuned initial condition
%

g = 1;
D = 1;
h = 0.8;

T = 10;
dt = 0.0001;
nt = round(T/dt)+1;
timey = linspace(0,T,nt);
x = zeros(nt,1);

% x(1) = 0;
x(1) = -log((0.75*exp(h)+0.25)/(0.25*exp(h)+0.75));

for j=1:nt-1,
   x(j+1)=x(j)+g*dt+sqrt(2*D*dt)*randn;
   if abs(x(j+1))>h, break, end
end

figure(1), plot(timey(1:j+1),x(1:j+1),'b','linewidth',3);
hold on, plot(timey(1:j+1000),h+0*timey(1:j+1000),'k','linewidth',2);
hold on, plot(timey(1:j+1000),-h+0*timey(1:j+1000),'k','linewidth',2);
set(gca,'xtick',[]); set(gca,'ytick',[]);
axis([0 timey(j+1000) -h-0.4 h+0.4])
box off

