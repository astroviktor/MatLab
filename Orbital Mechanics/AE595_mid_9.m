%% ORBITAL MECHANICS - MIDTERM
%% Ex. 9
clear all
close all
a=7000;             %semi major axis [km]
e=0.1;              %eccentricity
theta=[0:0.1:180];  %true anomaly range [degrees]
f=theta*pi/180;     %true anomaly in radians
mu=398600;          %gravitational parameter [km^3/s^2]

h=sqrt(a*mu*(1-e^2));

for i=1:length(theta)
    %position as a function of true anomaly
    r(i)=(h^2/mu)/(1+e*cos(f(i)));
    %velocity as a function of true anomaly (given energy equation)
    v(i)=sqrt(mu*((2/r(i))-(1/a)));
end

plot(theta,v,'LineWidth',1.5)
title('Orbital Velocity as a function of True Anomaly')
xlabel('\theta (°)')
ylabel('Velocity (km/s)')

    