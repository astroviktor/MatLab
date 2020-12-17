function [t]=time_elapsed_ellipsis(n,e,theta)

T=(2*pi)/n; %calculatin orbit period
mu=398600; %earth gravitational parameter

E=2*atan(sqrt((1+e)/(1-e))*tan(theta/2));    %computing eccentric anomaly
Me=E-e*sin(E);  %computin mean anomaly
t=(Me/2*pi)*T; %time elapsed since perigee at given theta

%Note: this function give the time elapsed from perigee just
%for elliptic orbits (0<e<1)

    
    
    
    