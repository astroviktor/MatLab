%% AE 504 - Compressible Flow
%% Homework 4
% Vittorio Baraldi

clear all; close all; clc

%% Ex.2
%% Part 1
%reading variables for Mach = 2.5 for the NACA 0012 airfoil
[x,y,Cp,P]=readvars('Mach25_NACA0012.txt');
%reading variables for Mach = 2.5 for the diamond-approximated airfoil
[x2,y2,Cp2,P2]=readvars('Mach25_Diamond.txt');
%mach number
M1=2.5;
%heat air coefficient
g=1.4;
%dynamic pressure for both results
q=g*M1^2*101325/2;
%angle of attack [rad]
alpha=3*pi/180;

%x and y coordinates for the quarter chord location. They are going to be
%needed for the pitching moment @ c/4 calculations later on
chordy=(y(223)+y(34))/2;
chordx=x(223);

%% NACA 0012
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
    if i>=129
        %if i>=129 then we are at the bottom surface, therefore the normal
        %force's contribution to the lift is positive
        L(i)=N(i)*cos(alpha);
    else
        %if i<129 then we are at the top surface, therefore the normal
        %force's contribution to the lift is positive
        L(i)=-N(i)*cos(alpha);
    end
    
    %Drag and Momentum
    if i>=129 && theta(i)<0
        %if i>=129 and dy/dx<0, then we are in the lower-left quadrant.
        %Knowing that we consequently know if the contribution of the
        %normal force to the drag is either positive or negative.
        %Same reasoning applies to the moment calculation
        D(i)=N(i)*sin(alpha);
        M(i)=-abs(L(i))*distanceL(i)+abs(D(i))*distanceD(i);
    elseif i<129 && theta(i)<0
        %if i<129 and dy/dx<0, then we are in the upper-right quadrant.
        D(i)=-N(i)*sin(alpha);
        M(i)=-abs(L(i))*distanceL(i)+abs(D(i))*distanceD(i);
    elseif i<129 && theta(i)>0
        %if i<129 and dy/dx>0, then we are in the upper-left quadrant.
        D(i)=N(i)*sin(alpha);
        M(i)=abs(L(i))*distanceL(i)-abs(D(i))*distanceD(i);
    elseif i>=129 && theta(i)>0
        %if i>=129 and dy/dx>0, then we are in the lower-right quadrant.
        D(i)=-N(i)*sin(alpha);
        M(i)=abs(L(i))*distanceL(i)-abs(D(i))*distanceD(i);
    end
end

%calculating the sumation of L, M and D arrays to obtain the total values
finalL=sum(L);
finalD=sum(D);
finalM=sum(M);
%computing lift, drag and moment coefficient knowing the dynamic pressure
%(S = 1)
cm=finalM/q;
cd=finalD/q;
cl=finalL/q;
%fixing theta, L, D, M arrays for plotting purposes
theta(1)=theta(2);L(1)=L(2);D(1)=D(2);M(1)=M(2);

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
%plotting profile vs Moment
figure(4)
plot(x*700,y*700,'-b','LineWidth',1.5)
hold on
plot(x*700,M,':r','LineWidth',2)
title('NACA 0012 vs Moment')
xlabel('^{x}/_{c} (/700)')
ylabel('Moment (Nm)')

%outputs
disp('-------------------------------------------------------------------')
fprintf('NACA 0012, real shape at alpha = 3°:\n')
fprintf('Lift coefficient: %1.4f\n',cl)
fprintf('Drag coefficient: %1.4f\n',cd)
fprintf('Pitching moment coefficient: %1.4f\n',cm)
fprintf('Drag: %1.4f N\n',finalD)
fprintf('Lift: %1.4f N\n',finalL)
disp('-------------------------------------------------------------------')

%% Diamond approximation

%In the following block of code the same principle as above is applied.
%Therefore comments are neglected in the following part (see above comments
%for references)
chordy2=(y2(354)+y2(52))/2;
chordx2=x2(354);
for i=2:length(x2)
    theta2(i)=atan((y2(i)-y2(i-1))/(x2(i)-x2(i-1)));
    l2(i)=sqrt((y2(i)-y2(i-1))^2+(x2(i)-x2(i-1))^2);
    ycoordinate2(i)=(y2(i)+y2(i-1))/2;
    xcoordinate2(i)=(x2(i)+x2(i-1))/2;
    distanceL2(i)=sqrt((xcoordinate2(i)-chordx2)^2+(ycoordinate2(i)-chordy2)^2);
    distanceD2(i)=abs(ycoordinate2(i)-chordy2);
    N2(i)=(P2(i)+P2(i-1))*l2(i)/2;
    if i>=206
        L2(i)=N2(i)*cos(alpha);
    else
        L2(i)=-N2(i)*cos(alpha);
    end
    if i>206 && theta2(i)<0
        D2(i)=N2(i)*sin(alpha);
        M2(i)=-abs(L2(i))*distanceL2(i)+abs(D2(i))*distanceD2(i);
    elseif i<=206 && theta2(i)<0
        D2(i)=-N2(i)*sin(alpha);
        M2(i)=-abs(L2(i))*distanceL2(i)+abs(D2(i))*distanceD2(i);
    elseif i<=206 && theta2(i)>0
        D2(i)=N2(i)*sin(alpha);
        M2(i)=abs(L2(i))*distanceL2(i)-abs(D2(i))*distanceD2(i);
    elseif i>206 && theta2(i)>0
        D2(i)=-N2(i)*sin(alpha);
        M2(i)=abs(L2(i))*distanceL2(i)-abs(D2(i))*distanceD2(i);
    end
end
finalL2=sum(L2);
finalD2=sum(D2);
finalM2=sum(M2);
cm2=finalM2/q;
cd2=finalD2/q;
cl2=finalL2/q;

theta2(1)=theta2(1);L2(1)=L2(1);D2(1)=D2(1);M2(1)=M2(1);

%PLOTS
figure(5)
plot(x2*50,y2*50,'-b','LineWidth',1.5)
hold on
plot(x2*50,theta2*180/pi,'-r','LineWidth',2)
title('NACA 0012 (diamond approximation) vs \theta')
xlabel('^{x}/_{c} (/50)')
ylabel('\theta (°)')
figure(6)
plot(x2*1000,y2*1000,'-b','LineWidth',1.5)
hold on
plot(x2*1000,L2,'-k','LineWidth',2)
title('NACA 0012 (diamond approximation) vs Lift')
xlabel('^{x}/_{c} (/1000)')
ylabel('Lift (N)')
figure(7)
plot(x2*100,y2*100,'-b','LineWidth',1.5)
hold on
plot(x2*100,D2,'-g','LineWidth',2)
title('NACA 0012 (diamond approximation) vs Drag')
xlabel('^{x}/_{c} (/100)')
ylabel('Drag (N)')
figure(8)
plot(x2*700,y2*700,'-b','LineWidth',1.5)
hold on
plot(x2*700,M2,':r','LineWidth',2)
title('NACA 0012 (diamond approximation) vs Moment')
xlabel('^{x}/_{c} (/700)')
ylabel('Moment (Nm)')

%outputs
disp('-------------------------------------------------------------------')
fprintf('NACA 0012, diamond approximation at alpha = 3°:\n')
fprintf('Lift coefficient: %1.4f\n',cl2)
fprintf('Drag coefficient: %1.4f\n',cd2)
fprintf('Pitching moment coefficient: %1.4f\n',cm2)
fprintf('Drag: %1.4f N\n',finalD2)
fprintf('Lift: %1.4f N\n',finalL2)
disp('-------------------------------------------------------------------')
