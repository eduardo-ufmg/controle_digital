function [x,pH,xc,pHc] = simrk_pH(x0,u1,u2,h,t,par,Kas,Ts)

xc = zeros(length(x0),Ts/h);
pHc = zeros(1,Ts/h);
xc(:,1) = x0;
p = roots([1 (Kas(1) - xc(2,1)) (Kas(1)*Kas(2)-Kas(1)*xc(2,1)-Kas(3)-Kas(1)*xc(3,1)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*xc(2,1)+2*Kas(1)*Kas(2)*xc(3,1)) -Kas(1)*Kas(2)*Kas(3)]);
ch = max(real(p));
pHc(1) = -log10(ch);

for i=2:Ts/h
    xc(:,i) = rkpH(xc(:,i-1),u1,u2,h,t,par);
    t = t+h;
    

    p = roots([1 (Kas(1) - xc(2,i)) (Kas(1)*Kas(2)-Kas(1)*xc(2,i)-Kas(3)-Kas(1)*xc(3,i)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*xc(2,i)+2*Kas(1)*Kas(2)*xc(3,i)) -Kas(1)*Kas(2)*Kas(3)]);
    ch = max(real(p));
    pHc(i) = -log10(ch); 
end

x = xc(:,i);

pH = pHc(i);

