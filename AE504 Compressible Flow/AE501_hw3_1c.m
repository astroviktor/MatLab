%% AE 504 Compressible Flow
%% Vittorio Baraldi
% Homework 3

clear all
close all
clc

%% Ex. 1
% Part c)

%initial data
Pres=3.79;          %reservoir pressure [atm]
Pe=1;               %exit pressure [atm]
At=0.5;             %throat area [m^2]
Ae=3.5;             %exit area [m^2]
Me=0.31;            %exit mach no.
g=1.4;              %heat air coefficient

exitratio=Ae/At;        

syms asol msol
A=2;                %initial shock position guess
dA=.001;            %increment
options = optimoptions('fsolve','Display','none');
for i=1:1000
    %area-mach no. function
    fM=@(x) ((((2/(g+1))*(1+0.2*x^2))^6)/x^2)-(A/At)^2;
    %solving for Mach no. before shock
    M1=fsolve(fM,2,options);
    %mach no. after shock
    M2=sqrt(((g-1)*M1^2+2)/(2*g*M1^2-(g-1)));
    %total pressure ratio across shock
    Pratio=((((g+1)*M1^2)/((g-1)*M1^2+2))^(g/(g-1)))*...
        (((g+1)/(2*g*M1^2-(g-1)))^(1/(g-1)));
    %Area ratio based on M2
    fA=@(x) ((((2/(g+1))*(1+0.2*M2^2))^6)/M2^2)-x^2;
    Aratio=fsolve(fA,1.1,options);
    %finding the new area ratio between exit area and the sonic area
    %condition after shock
    Ae_A2t=exitratio*(At/A)*Aratio;
    %solving for the exit Mach no. given the area ratio above
    fM2=@(x) ((((2/(g+1))*(1+0.2*x^2))^6)/x^2)-Ae_A2t^2;
    Me=fsolve(fM2,0.9,options);
    %with isentropic relations, find the conditions at the exit
    [Pexit,~,~,~,~,~]=isentropicflow(Me,g);
    %Pfinal is the exit pressure, given the initial guess
    Pfinal=Pexit*Pratio*Pres;
    %difference between computed exit pressure and actual exit pressure
    diff=abs(Pfinal-Pe);
    %checking on the difference tolerance (1e-3)
    if diff<1e-3
        sol=A;
        break
    else
        %if condition is not satisfied, guess is incremented
        A=A+dA;
    end
end

fprintf('Shock occurs at A = %1.2f meters squared', sol)
    
    
    
    
    