%% ORBITAL MECHANICS - MIDTERM
%% Ex. 13
clear all
close all
%from the equation given we can obtain the following parameters for the
%orbit

e=0.18;                     %eccentricity
p=15000;                    %parameter [km]
a=p/(1-e^2);                %semi major axis [km]
mu=398600;                  %gravitational parameter [km^3/s^2]
h=sqrt(p*mu);
%% Part a)

           
f=140*pi/180;
n=sqrt(mu/a^3);             %mean orbit rate [rad/s]
T=pi/n;                     %semi-period of the orbit [sec]

E = 2*atan( sqrt( (1-e)/(1+e) ) * tan(f/2) );
M = E - e*sin(E);
time = M/n;
t=T+(T-time);
fprintf('a) Time from perigee: %.0f seconds\n',t)

%% Part b)

t2=1000;                    %time [sec]

M2=n*t2; %mean anomaly
Ef=@(x) x-e*sin(x)-M2; %equation for Eccentric anomaly
dEf=@(x) 1-e*cos(x);   %derivative of the Eccentric anomaly eqn

E2=NewtonRhapsonSolver(Ef,dEf,M2,1e-10);

truean=2*atan( sqrt((1+e)/(1-e))*tan(E2/2) );

fprintf('b) True anomaly after 1000 seconds: %.0f°\n',truean*180/pi)
