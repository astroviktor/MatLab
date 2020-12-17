%% AE 504
% Homework 2

clear all
close all
clc

%% Ex. 1
% Part d)

M1=2.2;
g=1.4;

syms msol betasol thetasol

theta=15;

for i=1:500
    %wedge angle in radians
    thetarad=theta*pi/180;
    %wave angle from the theta-beta-m eqn
    beta1=abs(double(vpasolve(tan(thetarad)/(2*cot(betasol)) == ...
                ((M1^2*sin(betasol)^2)-1)/(M1^2*(g+cos(2*betasol))+...
                2),betasol,0.6071)));
             %normal mach number before shock
    M1n=M1*sin(beta1);
    %normal mach number after shock
    M2n=sqrt((M1n^2+(2/(g-1)))/(((2*g)/(g-1))*M1n^2-1));
    %mach number after shock
    M2=M2n/sin(beta1-thetarad);
    %mach angle for M2
    mu=asin(1/M2)*180/pi;
    %array of beta angles (from mu to pi/2)
    betaarray=[mu:0.5:90]*pi/180;
    for j=1:length(betaarray)
    theta2(j)=abs(double(vpasolve(tan(thetasol)/(2*cot(betaarray(j))) == ...
                ((M2^2*sin(betaarray(j))^2)-1)/(M2^2*(g+cos(2*betaarray(j)))+...
                2),thetasol,0.6071)));
    end
    thetamax=max(theta2);
    diff=thetamax-thetarad;
    %checking condition
    if diff<1e-3
        sol=M1;
        break
    else
        M1=M1-0.01;
        
    end
end

fprintf("The Mach reflection occurs at Mach = %1.2f",sol)
    
    