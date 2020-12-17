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
%Cl and Cd for NACA 0012
dyup=(y(129:257)-y(128:256))/(x(129:257)-x(128:256));
dylo=[];
cpup=[Cp(129:257,1)];
cplo=[];
for i=130:-1:2
    dylo=[dylo;(y(i,1)-y(i-1,1))/(x(i,1)-x(i-1,1))];
end
for i=129:-1:1
    cplo=[cplo;Cp(i,1)];
end
range=[x(1:129,1)];
cn=trapz(range,cpup-cplo);
ca0=trapz(range,((cplo.*dylo)-(cpup.*dyup)));
ca=trapz(range,ca0);
cl=cn*cos(alpha)-ca*sin(alpha);
cd=cn*sin(alpha)+ca*cos(alpha);
L=cl*q;
D=cd*q;
c=max(x)-min(x);
t=max(y)-min(y);
M=L*c/4+D*t/2;
cm=M/q/c;

%Cl and Cd for diamond approximation
dyup2=[];
for i=205:1:406
dyup2=[dyup2;(y2(i)-y2(i-1))/(x2(i)-x2(i-1))];
end
dyup2=[dyup2;(y2(406)-y2(1))/(x2(406)-x2(1))];
dylo2=[];
cpup2=[Cp2(205:406,1);Cp2(1,1)];
cplo2=[];
for i=204:-1:2
    dylo2=[dylo2;(y2(i,1)-y2(i-1,1))/(x2(i,1)-x2(i-1,1))];
end
for i=204:-1:2
    cplo2=[cplo2;Cp2(i,1)];
end
range2=[x2(1:203,1)];
cn2=trapz(range2,cpup2-cplo2);
ca2=trapz(range2,((cplo2.*dylo2)-(cpup2.*dyup2)));
cl2=cn2*cos(alpha)-ca2*sin(alpha);
cd2=cn2*sin(alpha)+ca2*cos(alpha);
L2=cl2*q;
D2=cd2*q;
c2=max(x2)-min(x2);
t2=max(y2)-min(y2);
M2=L2*c2/4+D2*t2/2;
cm2=M2/q/c2;
