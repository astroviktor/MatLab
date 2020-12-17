function [theta]=trueanomaly_ellipsis(time,e,n)

T=2*pi/n; %calculating orbit period
mu=398600; %earth gravitational parameter

Me=2*pi*(time/T);   %calculating mean anomaly through inverse equations
f=@(E) E-e*sin(E)-Me; %defining function for eccentric anomaly to be calculated
                       %through numerical methods
j=@(E) 1-e*cos(E);  %in order to use newton raphson method we need the
                    %derivative of the funcion for E, which is j(E)
toll=10^(-5);       %this represents the precision we want to obtain
                    %with the numerical method: it will not stop until
                    %this precision
maxit=500;          %max number of iteration allowed
x0=Me;              %value to start with: picked mean anomaly
E=newton_raphson(f,j,x0,toll,maxit)

theta=2*atan((tan(E/2))/(sqrt((1+e)/1-e)));

