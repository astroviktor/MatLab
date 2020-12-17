%Problem 3
%Homework 8
%AEM 4301 Orbital Mechanics

rm=2.28e8; %mars to sun distance (km)
re=1.496e8; %earth to sun distance (km)
mus=1.327e9; %sun gravitational parameter

at=0.5*(rm+re); %semimayor axis of the transfer orbit
nt=sqrt(mus/(at^3));
t12=pi/nt;
t=t12/2; %final transfer time

[v1,v2]=LambertSolver(re,rm,t,mus);



