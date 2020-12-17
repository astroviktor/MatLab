%ORBITAL MECHANICS 
%HOMEWORK 2
%Exercise 1

%earth gravitational parameter
mu=398600; %km^3/s^2

%CIRCULAR ORBIT
%initial conditions

r0=[39470;15939;9203]; %km
v0=[-1.279; 2.375;1.371];%km/s
v=norm(v0);
r=norm(r0);
T=2*pi*sqrt(r^3/mu);%orbit time

tspan=[0 T]; %sec
y0=[r0;v0];

%solving the differential equation of motion
    [t,y]=ode45(@orbitfun,tspan,y0);
    
    PlotEarth
    plot3(y(:,1),y(:,2),y(:,3))
    
%ELLIPTICAL ORBIT
    r0=[	-3859; 18952; 10942];%km
    v0=[-5.16;-0.263;-0.152];%km/s
    v=norm(v0);
    r=norm(r0);
    
    h=cross(r0,v0); %angular momentum vector and magnitude
    hn=norm(h); 
    e=0.5; %eccentricity
    a=(hn^2/mu)*(1/(1-e^2));  %semi-mayor axis
    T=2*pi*sqrt(a^3/mu);     %orbit time
    
    y0=[r0;v0];
    tspan=[0 T];
    [t,y]=ode45(@orbitfun,tspan,y0);
    
    figure(2)
    PlotEarth
    plot3(y(:,1),y(:,2),y(:,3))
    
%HYPERBOLIC ORBIT
    r0=[13015;986;569];
    v0=[-0.73;7.227;4.172];
    v=norm(v0);
    r=norm(r0);
       
    y0=[r0;v0]; 
    tspan=[0 T/3];
    [t,y]=ode45(@orbitfun,tspan,y0);
    
    figure(3)
    PlotEarth
    plot3(y(:,1),y(:,2),y(:,3))
    hold on
    plot3(y(:,1),-y(:,2),-y(:,3))






%--comments
%{
In order to know which conic section we have based on our initial
conditions we need to analyze different properties of the orbit such as
eccentricity (e), semi mayor axis (a), energy(E) and how velocity is
related to position. Each conic section has specific ranges/equations to
define such parameters and if we put them all together we can find a proper
guess for our initial position and velocity vectors.

In this problem, to make it easier and more accurate, an orbit online
simulator which will give a position and velocity vectors based on the
orbital elements. In order to obtain the requested conic sections, the
eccentricity value was changed. Credits http://en.homasim.com/orbitsimulation.php
%}

