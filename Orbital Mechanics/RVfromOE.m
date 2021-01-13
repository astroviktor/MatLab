function [r,v]=RVfromOE(a,e,i,omega,w,theta,mu)

h=sqrt(mu*abs(a*(1-e^2)));

rx = h^2/mu*(cos(theta)/(1-e*cos(theta));
ry= h^2/mu*(cos(theta)/(1-e*cos(theta));
vx=mu/h*e*sin(theta);
vy=mu/h*(e+cos(theta));
%r and v in perifocal frame
r=[rx,ry,0];
v=[vx,vy,0];

Q=perifocaltogeo(w,omega,i);
r=Q*r';
v=Q*v';




