%% AE 504 Compressible Flow
% Homework 1

clear all
close all

%% Ex.1
%Part a)
syms m1sol m2sol
%% Inviscid flow at M=0.8
T_inf=300;              %free-stream temperature [K]
P_inf=101325;           %free-stream pressure [Pa]
P1=P_inf;
M_inf=0.8;              %free-stream Mach no.
g=1.4;                  %heat air coefficent
R=287;
Cp=1.006;
a=sqrt(g*R*T_inf);
V=(M_inf*a);

%total temperature (free-stream condition)
T0_inf=Cp*T_inf+(V^2/2/1000);
%sonic temperature (free-stream condition)
Ts_inf=T0_inf*0.833;
%total pressure (free-stream condition)
P0_inf=P_inf*(T0_inf/T_inf)^(g/(g-1));
%sonic pressure (free-stream condition)
Ps_inf=P0_inf*0.528;

%from tecplot 
T1=275;
%from the aerodynamic calculator
M1=sqrt((2/(g-1))*(T0_inf/T1-1));
%after the shock
M2=sqrt(((g-1)*M1^2+2)/(2*g*M1^2-(g-1)));
P2=((g*M1^2*2-(g-1))/(g+1))*P1;
P01=P1/0.48379196;
P02=(((((g+1)*M1^2)/((g-1)*M1^2+2))^(g/(g-1)))*(((g+1)/...
    (2*g*M1^2-(g-1)))^(1/(g-1))))*P01;
T2=(((2*g*M1^2-(g-1))*((g-1)*M1^2+2))/((g+1)^2*M1^2))*T1;
%entropy change
dS=-R*log(P02/P01);
%theoretical mach numbers
Ms1=double(vpasolve(M_inf^2==2/(((g+1)/m1sol^2)-(g-1)),1.2));
Ms2=double(vpasolve(M2^2==2/(((g+1)/m2sol^2)-(g-1)),0.6));
Ms1=Ms1(real(Ms1)>0);
Ms2=Ms2(real(Ms2)>0);

fprintf("Conditions at point 2 for inviscid transonic flow:\n");
fprintf("Total pressure (P02) = %1.3f Pascal\n",P02);
fprintf("Static temperature (T2) = %1.3f K\n",T2);
fprintf("Mach number (M2) = %1.3f\n",M2);
fprintf("Theoretical Mach number 1 = %1.3f\n",Ms1);
fprintf("Theoretical Mach number 2 = %1.3f\n",Ms2);
fprintf("Entropy change (delta-S) = %1.3f\n\n",dS);

%% Viscous flow at M=0.8

%from tecplot 
T1=260;
%from the aerodynamic calculator
M1=sqrt((2/(g-1))*(T0_inf/T1-1));
%after the shock
M2=sqrt(((g-1)*M1^2+2)/(2*g*M1^2-(g-1)));
P2=((g*M1^2*2-(g-1))/(g+1))*P1;
P01=P1/0.48379196;
P02=(((((g+1)*M1^2)/((g-1)*M1^2+2))^(g/(g-1)))*(((g+1)/...
    (2*g*M1^2-(g-1)))^(1/(g-1))))*P01;
T2=(((2*g*M1^2-(g-1))*((g-1)*M1^2+2))/((g+1)^2*M1^2))*T1;
%entropy change
dS=-R*log(P02/P01);
%theoretical mach numbers
Ms1=double(vpasolve(M1^2==2/(((g+1)/m1sol^2)-(g-1)),1.2));
Ms2=double(vpasolve(M2^2==2/(((g+1)/m2sol^2)-(g-1)),0.6));
Ms1=Ms1(real(Ms1)>0);
Ms2=Ms2(real(Ms2)>0);

fprintf("Conditions at point 2 for viscous transonic flow:\n");
fprintf("Total pressure (P02) = %1.3f Pascal\n",P02);
fprintf("Static temperature (T2) = %1.3f K\n",T2);
fprintf("Mach number (M2) = %1.3f\n",M2);
fprintf("Theoretical Mach number 1 = %1.3f\n",Ms1);
fprintf("Theoretical Mach number 2 = %1.3f\n",Ms2);
fprintf("Entropy change (delta-S) = %1.3f\n\n",dS);

%% Inviscid flow at M=2
T_inf=300;              %free-stream temperature [K]
P_inf=101325;           %free-stream pressure [Pa]
P1=P_inf;
M_inf=2;              %free-stream Mach no.
g=1.4;                  %heat air coefficent
R=287;
Cp=1.006;
a=sqrt(g*R*T_inf);
V=(M_inf*a);

%from aerodynamic calculator
P01=P1/0.12780452;
P02=P1/0.17729110;
%from normal shock relations (assuming M_inf=M1)
M2=sqrt(((g-1)*M_inf^2+2)/(2*g*M_inf^2-(g-1)));
%from normal shock relations (assuming T_inf=T1)
T2=(((2*g*M_inf^2-(g-1))*((g-1)*M_inf^2+2))/((g+1)^2*M_inf^2))*T_inf;
dS=-R*log(P02/P01);
%theoretical mach numbers from eqn. 3.37
Ms1=double(vpasolve(M_inf^2==2/(((g+1)/m1sol^2)-(g-1)),1.5));
Ms2=double(vpasolve(M2^2==2/(((g+1)/m2sol^2)-(g-1)),0.5));
Ms1=Ms1(real(Ms1)>0);
Ms2=Ms2(real(Ms2)>0);

fprintf("Conditions at point 2 for inviscid supersonic flow:\n");
fprintf("Total pressure (P02) = %1.3f Pascal\n",P02);
fprintf("Static temperature (T2) = %1.3f K\n",T2);
fprintf("Mach number (M2) = %1.3f\n",M2);
fprintf("Theoretical Mach number 1 = %1.3f\n",Ms1);
fprintf("Theoretical Mach number 2 = %1.3f\n",Ms2);
fprintf("Entropy change (delta-S) = %1.3f\n\n",dS);


