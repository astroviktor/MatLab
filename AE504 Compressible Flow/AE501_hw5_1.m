%% AE 504 - Compressible Flow
%% Homework 5
clc;clear all;close all
%% Ex. 1
%initial conditions (L and R)
rhol=1;pl=1;ul=0;
rhor=0.15;pr=0.15;ur=0;g=1.4;
%time in seconds
t=0.2;
%starting location (L) in meters
x0=0;

%speed of sound (non-dimensional) (L)
aL=sqrt(g*pl/rhol);
%A and B are parameters used to simplify the expressions
A=(g+1)/(g-1);
B=2*g/(g+1);
%compatibility function
Ms=@(x) x-1/x-aL*A*(1-((pr/pl)*(B*x^2-1/A))^((g-1)/2/g));
Msol=fzero(Ms,1.2);
%pressure, density and velocity at point 1 and 2 (i.e. left and right of
%the contact region
p1=(B*Msol^2-1/A)*pr;
p2=p1;
rho1=rhor/(2/(g+1)/Msol^2+1/A);
u1=2/(g+1)*(Msol-1/Msol);
u2=u1;
a2=(p1/pl)^((g-1)/2/g)*aL;
rho2=(p2/pl)^(1/g)*rhol;
%location 1 and 2
x1=x0-aL*t;
x2=x0+(u2-a2)*t;
%x1 and x2 are the points which delimit the expansion region.
%In between that expansion:
x=linspace(x1,x2,30);
for i=1:length(x)
    u(i)=2/(g+1)*(aL+(x(i)-x0)/t);
    a(i)=aL-(g-1)*u(i)/2;
    p(i)=pl*(a(i)/aL)^(2*g/(g-1));
    rho(i)=(p(i)/pl)^(1/g)*rhol;
end
%x3 and x4
x3=x0+u2*t;
x4=x0+Msol*t;
%concatenating arrays for plot
xarray=[-0.5;x1;x';x2;x3;x3;x4;x4;0.5];
rhoarray=[rhol;rhol;rho';rho2;rho2;rho1;rho1;rhor;rhor];
parray=[pl;pl;p';p2;p2;p1;p1;pr;pr];
uarray=[ul;ul;u';u2;u2;u1;u1;ur;ur];
%outputs and plots
figure(1)
plot(xarray,rhoarray,'-r','LineWidth',2)
title('Plot Density vs Position')
grid on
xlabel('x(m)')
ylabel('Density(kg/m^3)')
figure(2)
plot(xarray,parray,'-b','LineWidth',2)
title('Plot Pressure vs Position')
grid on
xlabel('x(m)')
ylabel('Pressure (Pa)')
figure(3)
plot(xarray,uarray,'-g','LineWidth',2)
title('Plot Velocity vs Position')
grid on
xlabel('x(m)')
ylabel('Velocity (m/s)')
