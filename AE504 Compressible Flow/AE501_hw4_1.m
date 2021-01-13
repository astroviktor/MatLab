%% AE 504 - Compressible Flow
%% Homework 4
% Vittorio Baraldi

clear all; close all; clc

%% Ex. 1
%% Part 1
%file reading
[x,y,Cp,P]=readvars('NACA_0012_a3_inc.txt');
[x3,y3,Cp3,P3]=readvars('NACA_0012_a3_Mach_0.5.txt');
%mach number
M1=0.5;
%Prandtl-Glauert approximation
Cp_compressible=Cp/sqrt(1-M1^2);

%plot
figure(1)
plot(x3,Cp3,'-r','LineWidth',1.5)
hold on
plot(x,Cp_compressible,'-b','LineWidth',1.5)
legend('Pressure coefficient for compressible flow','Pressure coefficient approximation from inc. data for compressible flow')
xlabel('x')
ylabel('C_p')
title('Pressure coefficient C_p vs Chord location')
hold off

%% Part 2
syms msol
%Critical pressure coefficient (minimum Cp)
Cp_cr=min(Cp);
g=1.4;          %heat air coefficient
%critical mach no. 
M_cr=double(vpasolve(Cp_cr==(2/g/msol^2)*((((1+((g-1)/2)*msol^2)/...
    (1+(g-1)/2))^(g/(g-1)))-1),msol,0.5));

%outputs
disp('Critical Mach number:')
disp(M_cr)

%% Part 3
%file reading
[x2,y2,Cp2,P2]=readvars('NACA_0012_a3_Mach_0.8.txt');
%mach number for drag-M divergence
M2=linspace(0,1,200);
%dynamic pressure assuming M=0.1 for this example, since we are analyzing
%an incompressible airfoil
q=g*0.1^2*101325/2;
%angle of attack in radians
alpha=3*pi/180;

%x and y coordinates for the quarter chord location. They are going to be
%needed for the pitching moment @ c/4 calculations later on
chordy=(y(135)+y(67))/2;
chordx=x(67);
for i=2:length(x)
    %dy/dx (or theta angle)
    theta(i)=atan((y(i)-y(i-1))/(x(i)-x(i-1)));
    %panel lenght
    l(i)=sqrt((y(i)-y(i-1))^2+(x(i)-x(i-1))^2);
    %x and y coordinate of panel midpoint (for pitching moment)
    ycoordinate(i)=(y(i)+y(i-1))/2;
    xcoordinate(i)=(x(i)+x(i-1))/2;
    %distance of the lift force vector from the quarter chord
    distanceL(i)=sqrt((xcoordinate(i)-chordx)^2+(ycoordinate(i)-chordy)^2);
    %distance of the drag vector from the chord line
    distanceD(i)=abs(ycoordinate(i)-chordy);
    %normal force to the panel
    N(i)=(P(i)+P(i-1))*l(i)/2;
    
    %Lift
    if i>=101
        %if i>=129 then we are at the bottom surface, therefore the normal
        %force's contribution to the lift is positive
        L(i)=N(i)*cos(alpha);
    else
        %if i<129 then we are at the top surface, therefore the normal
        %force's contribution to the lift is positive
        L(i)=-N(i)*cos(alpha);
    end
    
    %Drag and Momentum
    if i>=101 && theta(i)<0
        %if i>=101 and dy/dx<0, then we are in the lower-left quadrant.
        %Knowing that we consequently know if the contribution of the
        %normal force to the drag is either positive or negative.
        %Same reasoning applies to the moment calculation
        D(i)=-N(i)*sin(alpha);
    elseif i<101 && theta(i)<0
        %if i<101 and dy/dx<0, then we are in the upper-right quadrant.
        D(i)=N(i)*sin(alpha);
    elseif i<101 && theta(i)>0
        %if i<101 and dy/dx>0, then we are in the upper-left quadrant.
        D(i)=-N(i)*sin(alpha);
    elseif i>=101 && theta(i)>0
        %if i>=101 and dy/dx>0, then we are in the lower-right quadrant.
        D(i)=N(i)*sin(alpha);
    end
end

%calculating the sumation of L, M and D arrays to obtain the total values
finalL=sum(L);
finalD=sum(D);
%computing lift, drag and moment coefficient knowing the dynamic pressure
%(S = 1)
cd=finalD/q;
cl=finalL/q;
%fixing theta, L, D, M arrays for plotting purposes
theta(1)=theta(2);L(1)=L(2);D(1)=D(2);

%PLOTS
%plotting the profile vs dy/dx
figure(1)
plot(x*100,y*100,'-b','LineWidth',1.5)
hold on
plot(x*100,theta*180/pi,'-r','LineWidth',2)
title('NACA 0012 vs \theta')
xlabel('^{x}/_{c} (/100)')
ylabel('\theta (°)')
%plotting profile vs Lift
figure(2)
plot(x*6000,y*6000,'-b','LineWidth',1.5)
hold on
plot(x*6000,L,'-k','LineWidth',2)
title('NACA 0012 vs Lift')
xlabel('^{x}/_{c} (/6000)')
ylabel('Lift (N)')
%plotting profile vs Drag
figure(3)
plot(x*100,y*100,'-b','LineWidth',1.5)
hold on
plot(x*100,D,'-g','LineWidth',2)
title('NACA 0012 vs Drag')
xlabel('^{x}/_{c} (/100)')
ylabel('Drag (N)')

%outputs
disp('-------------------------------------------------------------------')
fprintf('NACA 0012 incompressible coefficients = 3°:\n')
fprintf('Lift coefficient: %1.4f\n',cl)
fprintf('Drag coefficient: %1.4f\n',cd)
fprintf('Pitching moment coefficient: %1.4f\n',cm)
fprintf('Drag: %1.4f N\n',finalD)
fprintf('Lift: %1.4f N\n',finalL)
disp('-------------------------------------------------------------------')

%last part
M2=linspace(0,1,500);
for i=1:length(M2)
    Cd2(i)=(2*pi*alpha^2)/sqrt(1-M2(i)^2);
end

%plot 
figure(4)
plot(M2,Cd2,'-g','LineWidth',2)
xlabel('Mach number')
ylabel('C_d')
title('Mach number-drag divergence at \alpha = 3°')    

