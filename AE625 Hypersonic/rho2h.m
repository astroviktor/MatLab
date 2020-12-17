function [altitude]=rho2h(rho)  %the name of the function is rho2h
    
a0=-0.003567; %gradient input °R/ft
%gradient is the temperature increasing/decreasing per feet
g=32.2; %ft/sec^2 gravity acceleration
R=1716; %gas costant
rho0=2.37*(1e-3);

T0=518.67;%if we are in troposhpere or stratosphere calculations will change
esp=-(1+g/R*a0); %precalculating exponent to avoid long expressions
T=T0*((rho/rho0)^(1/esp));  %inverse formula from density
altitude=(T-T0)/a0;

    

