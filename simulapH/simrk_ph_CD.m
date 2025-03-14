clearvars; close all; clc

setup_pH;

h = 10;

t0 = h;
tm = 50;
tf = 60*tm;
tr = 60*15;
td = 60*10;
t = t0:h:tf;

Ts = 40;
T = t(1:Ts/h:end);

K = 2;

Q1 = 3*ones(1,length(T)); 

d1 = zeros(1,length(T)); 

Q3 = 2*ones(1,length(T)); 

r = 5*ones(1,length(T));

xc = zeros(length(x0),length(t)); 
pHc = zeros(1,length(t));

x = zeros(length(x0),length(T));
pH = zeros(1,length(T));

e = zeros(1,length(T));

m = zeros(1,length(T));

x(:,1) = x0;

for k = 2:floor(tr/Ts)+1
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),Q1(k),Q3(k),h,t(kc),par,Kas,Ts);

end

ini = k;

u1 = Q1;
u2 = Q3;

m = zeros(size(Q3));

r(ini:end) = 5.5*ones(1,length(T)-ini+1); 

d1(ini+floor(td/Ts):end) = 0.3*ones(1,length(T)-ini-floor(td/Ts)+1); 

for k=ini:length(T)

    e(k-1) = r(k-1)-pH(k-1);
      

    m(k-1) = K*e(k-1);
    u2(k-1) = Q3(k-1)+m(k-1);
    

    u1(k-1) = Q1(k-1) - d1(k-1);
    

    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),u1(k-1),u2(k-1),h,t(kc),par,Kas,Ts);
end

pplot = 1;
if pplot == 1

   figure(1)
   yyaxis left
   plot(t/60,pHc,'k',(T+Ts)/60,pH,'b*');
   set(gca,'FontSize',16)
   ylabel('pH')
   yyaxis right
   stairs((T+Ts)/60,u2,'r');
   set(gca,'FontSize',16)
   xlabel('t (min)')
   ylabel('vazao de base (mL/seg)')
   axis([0 (t(end)+h)/60 0 12])
end

