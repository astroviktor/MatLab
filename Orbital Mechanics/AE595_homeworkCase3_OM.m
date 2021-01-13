%CASE 3

r1 = [506.717068361234;-470.869269481925;138.096930726997];
r2 = [-10298.3068750911;6203.57304746482;-16637.5263693708];
dT = 1800.000000;

mu=398600.44;

[v1,v2,theta]=LambertSolver(r1,r2,dT,mu,'pro');
[v12,v22,theta2]=LambertSolver(r1,r2,dT,mu,'retro');
v1
v2
v12
v22
[a,e,i,omega,theta,w]=OEfromRV(r1,v1,mu);
[a2,e2,i2,omega2,theta2,w2]=OEfromRV(r1,v12,mu);
fprintf('Prograde Orbit:\nEccentricity: %1.2f\nInclination: %1.2f °\nArgument of Periapsis: %1.2f °\nLongitude of the Ascending node: %1.2f °\nTrue anomaly: %1.2f °\nSemi Mayor Axis: %4.2f km\n\n',e,i*180/pi,w*180/pi,omega*180/pi,theta*180/pi,a);
fprintf('Retrograde Orbit:\nEccentricity: %1.2f\nInclination: %1.2f °\nArgument of Periapsis: %1.2f °\nLongitude of the Ascending node: %1.2f °\nTrue anomaly: %1.2f °\nSemi Mayor Axis: %4.2f km\n',e2,i2*180/pi,w2*180/pi,omega2*180/pi,theta2*180/pi,a2);

if( e<1 )
  vtheta=linspace(0,2*pi,300);
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

