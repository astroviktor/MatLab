%CASE 2

r1 = [-5192.64707332683;-2730.67520032545;-1312.84187347459];
r2 = [-1396.17652892469;-2398.88103277394;7822.986845429];
dt = 1288.841280;
mu=398600.44;

[v1,v2,theta]=LambertSolver(r1,r2,dt,mu,'pro');
[v12,v22,theta2]=LambertSolver(r1,r2,dt,mu,'retro');
[a,e,i,omega,theta,w]=OEfromRV(r1,v1,mu);
[a2,e2,i2,omega2,theta2,w2]=OEfromRV(r1,v12,mu);
v1
v2
v12
v22
fprintf('Prograde Orbit:\nEccentricity: %1.2f\nInclination: %1.2f °\nArgument of Periapsis: %1.2f °\nLongitude of the Ascending node: %1.2f °\nTrue anomaly: %1.2f °\nSemi Mayor Axis: %4.2f km\n\n',e,i*180/pi,w*180/pi,omega*180/pi,theta*180/pi,a);
fprintf('Retrograde Orbit:\nEccentricity: %1.2f\nInclination: %1.2f °\nArgument of Periapsis: %1.2f °\nLongitude of the Ascending node: %1.2f °\nTrue anomaly: %1.2f °\nSemi Mayor Axis: %4.2f km\n',e2,i2*180/pi,w2*180/pi,omega2*180/pi,theta2*180/pi,a2);

if( e<1 )
  vtheta=linspace(0,2*pi,500);
else
  thetaMax=acos(-1/e);
  thetaMin=-thetaMax;
  vtheta=linspace(thetaMin+0.01,thetaMax-0.01,500);
end

if( e2<1 )
  vtheta2=linspace(0,2*pi,500);
else
  thetaMax2=acos(-1/e);
  thetaMin2=-thetaMax2;
  vtheta2=linspace(thetaMin2+0.01,thetaMax2-0.01,500);
end

PlotTraj(a,e,w,i,omega,theta,vtheta);
fig=gcf;
PlotTraj(a2,e2,w2,i2,omega2,theta2,vtheta2,fig);

