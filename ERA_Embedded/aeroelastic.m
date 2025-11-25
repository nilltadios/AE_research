function dxdt= aeroelastic(t,x,gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,Us)
dxdt=zeros(size(x));

Wratio=0.2;  %frequency ration of the airfoil
phi1=0.165; phi2=0.335; eps1=0.0455; eps2=0.3; mu=100; ah=-1/2; xA=1/4; zetaXi=0; zetaA=0; rA=.7;  %  

c0=1+1/mu; c1=xA-ah/mu; c2=2/mu*(1-phi1-phi2); 
c3=(1+2*(.5-ah)*(1-phi1-phi2))/mu; c4=gammaX*(Wratio/Us)^2+2/mu*(eps1*phi1+eps2*phi2);c5=(Wratio/Us)^2*gammaXXX;
c6=2/mu*(1-phi1-phi2+(1/2-ah)*(eps1*phi1+eps2*phi2));
c7=2/mu*eps1*phi1*(1-eps1*(1/2-ah)); c8=2/mu*eps2*phi2*(1-eps2*(1/2-ah));
c9=-2/mu*eps1^2*phi1; c10=-2/mu*eps2^2*phi2;

d0=xA/rA^2-ah/(mu*rA^2); d1=1+(1+8*ah^2)/(8*mu*rA^2);
d2=(1-2*ah)/(2*mu*rA^2)-(1+2*ah)*(1-2*ah)*(1-phi1-phi2)/(2*mu*rA^2);
d3=gammaA/Us^2-(1+2*ah)*(1-phi1-phi2)/(mu*rA^2)-(1+2*ah)*(1-2*ah)*(phi1*eps1+phi2*eps2)/(2*mu*rA^2);
d4=gammaAAA/Us^2;
d5=-(1+2*ah)*(1-phi1-phi2)/(mu*rA^2);
d6=-(1+2*ah)*(phi1*eps1+phi2*eps2)/(mu*rA^2);
d7=-(1+2*ah)*phi1*eps1*(1-eps1*(1/2-ah))/(mu*rA^2);
d8=-(1+2*ah)*phi2*eps2*(1-eps2*(1/2-ah))/(mu*rA^2);
d9=(1+2*ah)*phi1*eps1^2/(mu*rA^2);
d10=(1+2*ah)*phi2*eps2^2/(mu*rA^2);


% P=c2*x(4)+c3*x(2)+c4*x(3)+c5*x(3)^3+c6*x(1)+c7*x(5)+c8*x(6)+c9*x(7)+c10*x(8);
% H=d2*x(2)+d3*x(1)+d4*x(1)^3+d5*x(4)+d6*x(3)+d7*x(5)+d8*x(6)+d9*x(7)+d10*x(8);

P=c2*x(4)+c3*x(2)+c4*x(3)+c5*x(3)^3+gammaXXXXX*(Wratio/Us)^2*x(3)^5+c6*x(1)+c7*x(5)+c8*x(6)+c9*x(7)+c10*x(8)-2/mu*((1/2-ah)*x(1,1)+x(3,1))*(phi1*eps1*exp(-eps1*t)+phi2*eps2*exp(-eps2*t));
H=d2*x(2)+d3*x(1)+d4*x(1)^3+gammaAAAAA*(1/Us^2)*x(1)^5+d5*x(4)+d6*x(3)+d7*x(5)+d8*x(6)+d9*x(7)+d10*x(8)+(1+2*ah)/(2*rA^2)*(2/mu*((1/2-ah)*x(1,1)+x(3,1))*(phi1*eps1*exp(-eps1*t)+phi2*eps2*exp(-eps2*t)));


dxdt(1)=x(2);
dxdt(2)=(c0*H-d0*P)/(d0*c1-c0*d1);
dxdt(3)=x(4);
dxdt(4)=(-c1*H+d1*P)/(d0*c1-c0*d1);
dxdt(5)=x(1)-eps1*x(5);
dxdt(6)=x(1)-eps2*x(6);
dxdt(7)=x(3)-eps1*x(7);
dxdt(8)=x(3)-eps2*x(8);
