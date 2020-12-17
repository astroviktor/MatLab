%% AE 504 Compressible Flow
% Homework 1
clear all
close all

%% Ex. 2
T_inf=200;                      %free-stream temperature [K]
P_inf=0.057*101325;             %free-stream pressure [Pa]
M1=4.197;                       %initial mach number
g=1.4;                          %heat air coefficient

%total pressure at point 1 (ratio from aerodynamic calculator)
P01=P_inf/0.00508182;
%total temperature at point 1
T01=T_inf/0.22109406;
%static pressure at point 1
P1=P01*((1+((g-1)/2)*M1^2)^(-g/(g-1)));
%static temperature at point 1
T1=T01*(1/(1+((g-1)/2)*M1^2));

fprintf("Conditions at point 1:\n");
fprintf("Total pressure = %1.3f Pascal\n",P01);
fprintf("Total Temperature = %1.3f K\n",T01);
fprintf("Static pressure = %1.3f Pascal\n",P1);
fprintf("Static temperature = %1.3f K\n",T1);
fprintf("Mach number = %1.3f\n\n",M1);

%conditions at point 2
%after the shock
M2=sqrt(((g-1)*M1^2+2)/(2*g*M1^2-(g-1)));
P2=((g*M1^2*2-(g-1))/(g+1))*P1;
T2=(((2*g*M1^2-(g-1))*((g-1)*M1^2+2))/((g+1)^2*M1^2))*T1;
P02=(((((g+1)*M1^2)/((g-1)*M1^2+2))^(g/(g-1)))*(((g+1)/...
    (2*g*M1^2-(g-1)))^(1/(g-1))))*P01;
T02=T01;

fprintf("Conditions at point 2:\n");
fprintf("Total pressure = %1.3f Pascal\n",P02);
fprintf("Total Temperature = %1.3f K\n",T02);
fprintf("Static pressure = %1.3f Pascal\n",P2);
fprintf("Static temperature = %1.3f K\n",T2);
fprintf("Mach number = %1.3f\n\n",M2);


%conditions at point 3 (sonic conditions)
%ratios obtained from aerodynamic calculator at M2
T3=T2/1.15720650;
P3=P2/1.66700932;
P03=P3/((2/(g+1))^(g/(g-1)));
T03=T3*(g+1)/2;
M3=1;
fprintf("Conditions at point 3:\n");
fprintf("Total pressure = %1.3f Pascal\n",P03);
fprintf("Total Temperature = %1.3f K\n",T03);
fprintf("Static pressure = %1.3f Pascal\n",P3);
fprintf("Static temperature = %1.3f K\n",T3);
fprintf("Mach number = %1.3f\n\n",M3);


%conditions at point 4
%ratios obtained from aerodynamic calculator at M4
M4=2.197;
P4=P3*0.17786298;
T4=T3*0.61057460;
P04=P4/0.09396177;
T04=T4/0.50881216;

fprintf("Conditions at point 4:\n");
fprintf("Total pressure = %1.3f Pascal\n",P04);
fprintf("Total Temperature = %1.3f K\n",T04);
fprintf("Static pressure = %1.3f Pascal\n",P4);
fprintf("Static temperature = %1.3f K\n",T4);
fprintf("Mach number = %1.3f\n\n",M4);

%conditions at point 5
%after the shock
M5=sqrt(((g-1)*M4^2+2)/(2*g*M4^2-(g-1)));
P5=((g*M4^2*2-(g-1))/(g+1))*P4;
T5=(((2*g*M4^2-(g-1))*((g-1)*M4^2+2))/((g+1)^2*M4^2))*T4;
P05=(((((g+1)*M4^2)/((g-1)*M4^2+2))^(g/(g-1)))*(((g+1)/...
    (2*g*M4^2-(g-1)))^(1/(g-1))))*P04;
T05=T04;

fprintf("Conditions at point 5:\n");
fprintf("Total pressure = %1.3f Pascal\n",P05);
fprintf("Total Temperature = %1.3f K\n",T05);
fprintf("Static pressure = %1.3f Pascal\n",P5);
fprintf("Static temperature = %1.3f K\n",T5);
fprintf("Mach number = %1.3f\n\n",M5);

%conditions at point 6 (combustor)
T6=895;
T06=T05;
tempratio=T6/T06;
%from the tables, knowing T6/T06 we can get the mach no. at point 6
M6=0.23144722;
P06=P05;
P6=0.96338897*P06;

fprintf("Conditions at point 6:\n");
fprintf("Total pressure = %1.3f Pascal\n",P06);
fprintf("Total Temperature = %1.3f K\n",T06);
fprintf("Static pressure = %1.3f Pascal\n",P6);
fprintf("Static temperature = %1.3f K\n",T6);
fprintf("Mach number = %1.3f\n\n",M6);

%plotting
Parray=[P2 P3 P4 P5 P6];
Tarray=[T2 T3 T4 T5 T6];
x=[-2 -1.5 0 1.5 4];
figure(1)
plot(x,Parray,'-r')
xlabel('x')
ylabel('Pressure (Pa)')
hold on
figure(2)
plot(x,Tarray,'-b')
xlabel('x')
ylabel('Temperature (K)')
ylim([0 1000])

%{
Legend
M#     -- Mach number at point #
P0#    -- Total pressure at point #
P#     -- Static pressure at point #
T#     -- Static temperature at point #
T0#    -- Total temperature at point #
%}

%% Ex. 3

f=0.005;            %friction coefficient
D=0.5;              %tube diameter [m]
L=20;               %tube length [m]
M=0.4;              %mach no.
P=101325;           %pressure [Pa]
T=300;              %temperature [K]

%from eqn 2.80
%from friction tables
k1=2.308;
P0=P/0.04047715;
k2=k1-(4*f*L/D);
Ts=T/1.163;
Ps=P/2.696;
P0s=P0/1.590;

M_2=0.45; %approximation from tables
T_2=Ts*1.153;
P_2=Ps*2.383;
P_02=P0s*1.45;

fprintf("Conditions at the exit are:\n");
fprintf("Total pressure (P02) = %1.3f Pascal\n",P_02);
fprintf("Temperature (T2) = %1.3f K\n",T_2);
fprintf("Static pressure (P2) = %1.3f Pascal\n",P_2);


