function xdot = dvpH(x,ux,uy,t,par)

xd(1) = (ux+par(1)+uy-(par(9)*(sqrt(x(1)-par(10)))))/par(3); 

xd(2) = ((ux*(par(11)-x(2)))+(par(1)*(par(12)-x(2)))+(uy*(par(13)-x(2))))/(par(4)); 

xd(3) = ((ux*(par(14)-x(3)))+(par(1)*(par(15)-x(3)))+(uy*(par(16)-x(3))))/(par(4));

xd(4) = -ux/par(5); 

xd(5) = -par(1)/par(6); 

xd(6) = -uy/par(7); 

xd(7) = par(2)/par(8); 

xdot=xd';
end
