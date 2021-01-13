%% AE504 
%% Vittorio Baraldi
% Homework 2
clear all
close all
clc

%% Ex. 3
% Part 2

Mu=5.3;                     %mach number upper zone
Ml=2.2549;                  %mach number lower zone
g=1.4;                      %heat air coefficient
Pu=4784.69;
Pl=282156.232;

%In the upper zone, a shock wave occurs, while in the lower one we have an
%expansion wave

dangle=0.001;
delta=-3;
syms betasol msol
for i=1:10000
    %shock wave
    
    %wedge angle
    thetau=(30-delta)*pi/180;
    %wave angle from the theta-beta-m eqn
    beta=abs(double(vpasolve(tan(thetau)/(2*cot(betasol)) == ...
                ((Mu^2*sin(betasol)^2)-1)/(Mu^2*(g+cos(2*betasol))+...
                2),betasol,0.6071)));
    %normal mach number before shock
    M1n=Mu*sin(beta);
    %normal mach number after shock
    M2n=sqrt((M1n^2+(2/(g-1)))/(((2*g)/(g-1))*M1n^2-1));
    %mach number after shock
    M2u=M2n/sin(beta-thetau);
    %pressure in the upper zone
    P4up=(1+((2*g)/(g+1))*(M1n^2-1))*Pu;
    P4up=P4up/101325;
    
    %expansion wave
    
    %turn angle
    thetal=15-delta;
    %prandtl-meyer function before and after expansion
    v1=sqrt((g+1)/(g-1))*atand(sqrt(((g-1)/(g+1))*(Ml^2-1)))-...
        atand(sqrt(Ml^2-1));
    v2=v1+thetal;
    %mach number after expansion
    M2l=double(vpasolve(v2==sqrt((g+1)/(g-1))*atand(sqrt(((g-1)/(g+1))*...
        (msol^2-1)))-atand(sqrt(msol^2-1)),msol,1.4));
    %pressure in the lower zone
    P4lo=Pl/(((1+((g-1)/2)*M2l^2)/(1+((g-1)/2)*Ml^2))^(g/(g-1)));
    P4lo=P4lo/101325;
    
    %checking on difference pressure
    %if it gets close enough to zero, we have a result
    diff=abs(P4lo-P4up);
    if diff<1e-3
        %solution output variables
        deltasol=delta;
        M4up=M2u;
        M4lo=M2l;
        break
    else
        %if the difference in pressure does not satisfy the condition, we
        %will increase delta
        delta=delta+dangle;
    end
end

%outputs
fprintf("Flow direction (delta angle) = %1.3f°\n",deltasol);
fprintf("Mach number (upper zone) = %1.3f\n",M4up);
fprintf("Mach number (lower zone) = %1.3f\n",M4lo);

    