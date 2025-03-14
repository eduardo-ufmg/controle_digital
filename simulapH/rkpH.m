function [x,xc] = rkpH(x0,ux,uy,h,t,par)

xd = dvpH(x0,ux,uy,t,par);
savex0 = x0;
phi = xd;
x0 = savex0+0.5*h*xd;

xd = dvpH(x0,ux,uy,t+0.5*h,par);
phi = phi+2*xd;
x0 = savex0+0.5*h*xd;

xd = dvpH(x0,ux,uy,t+0.5*h,par);
phi = phi+2*xd;
x0 = savex0+h*xd;

xd = dvpH(x0,ux,uy,t+h,par);
x = savex0+(phi+xd)*h/6;
end