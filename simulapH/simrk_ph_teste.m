clearvars; close all; clc

setup_pH;

h = 10;

t0 = h;
tm = 50;
tf = 60*tm;
t = t0:h:tf;

Ts = 40;
T = t(1:Ts/h:end);

Q1 = 3*ones(1,length(T)); 

Q3 = 2*ones(1,length(T)); 
u = zeros(1,length(T)); 

xc = zeros(length(x0),length(t)); 
pHc = zeros(1,length(t));

x = zeros(length(x0),length(T));
y1 = zeros(1,length(T));
y2 = y1;
y3 = y1;
y4 = y1;
y5 = y1;
pH = zeros(1,length(T));

x(:,1) = x0;

for k = 2:600/Ts+1
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),Q1(k),Q3(k),h,t(kc),par,Kas,Ts);

end

for k = 6:600/Ts+1
    y1(k) = 0.7770*y1(k-1) + 0.6177*u(k-1); 
    y2(k) = 0.8521*y2(k-1) + 0.6234*u(k-1) -0.4327*u(k-3) + 0.2190*u(k-4);
    y3(k) = 1.4042*y3(k-1) -0.4983*y3(k-2) +0.3569*u(k-1) -0.5219*u(k-4) +0.4256*u(k-5);
    y4(k) = 0.8888*y4(k-1) +0.5221*u(k-1) -0.4247*u(k-4) +0.2106*u(k-5);
    y5(k) = 1.3316*y5(k-1) -0.4182*y5(k-2) +0.2399*u(k-1);
end

ini = k;

u1 = Q1;
u2 = Q3;

du2 = zeros(size(Q3));

du2(ini:end) = 0.5*ones(1,length(T)-ini+1); 
u2 = Q3 + du2;

for k=ini:length(T)
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),u1(k),u2(k),h,t(kc),par,Kas,Ts);
end

u = du2;
for k = ini:length(T)
   y1(k) = 0.7770*y1(k-1) + 0.6177*u(k-1); 
    y2(k) = 0.8521*y2(k-1) + 0.6234*u(k-1) -0.4327*u(k-3) + 0.2190*u(k-4);
    y3(k) = 1.4042*y3(k-1) -0.4983*y3(k-2) +0.3569*u(k-1) -0.5219*u(k-4) +0.4256*u(k-5);
    y4(k) = 0.8888*y4(k-1) +0.5221*u(k-1) -0.4247*u(k-4) +0.2106*u(k-5);
    y5(k) = 1.3316*y5(k-1) -0.4182*y5(k-2) +0.2399*u(k-1);
end

pplot = 1;
if pplot == 1

   figure(1)
   yyaxis left
   plot(t/60,pHc,'k',(T+Ts)/60,pH,'b*');
   set(gca,'FontSize',16)
   ylabel('pH')
   yyaxis right
   plot((T+Ts)/60,u2,'r');
   set(gca,'FontSize',16)
   xlabel('t (min)')
   ylabel('vazao de base (mL/seg)')
   axis([0 (t(end)+h)/60 0 12])
end

figure(2)
plot((T+Ts)/60,pH,'b',T/60,y1+5,'k+',T/60,y2+5,'r+',T/60,y3+5,'g+',T/60,y4+5,'b+',T/60,y5+5,'c+');
axis([9 T(end)/60 4.9 6.5])
set(gca,'FontSize',16)
ylabel('pH')
xlabel('t (min)')

