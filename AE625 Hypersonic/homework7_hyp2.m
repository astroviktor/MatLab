%% EX. 4.2
clear all
close all
Q1=120000000;                   %heat of reaction (hydrogen) [J/kg]
Q2=50010000;                    %heat of reaction (methane) [J/kg]
f=0.0291;                       %fuel/air ratio
g0=9.81;                        %acceleration of gravity [m/s^2]
Ve1=[1600:100:3000];             %exit velocity [m/s]
M=5;                            %mach numbers
a=303.1;                        %speed of sound at 30,000 feet [m/s]
V0=M*a;                         %flight speed
%% HYDROGEN

for j=1:length(Ve1)
    %dimensionless parameter
     k1(j)=(f*Q1)/(V0^2/2);  
     %thermal efficiency
     etath1(j)=((Ve1(j)^2/2)-(V0^2/2))/(f*Q1);
     %overall efficiency
     eta01(j)=(2*(sqrt(etath1(j)*k1(j)+1)-1))/k1(j);
     %specific impulse
     I1(j)=(Q1*eta01(j))/(g0*V0);
end

%% METHANE
Ve2=linspace(1550,2250,length(Ve1)); 
for j=1:length(Ve1)
    k2(j)=(f*Q2)/(V0^2/2);  
    etath2(j)=((Ve2(j)^2/2)-(V0^2/2))/(f*Q2);
    eta02(j)=(2*(sqrt(etath2(j)*k2(j)+1)-1))/k2(j);
    I2(j)=(Q2*eta02(j))/(g0*V0);
end

%% PLOTS
figure(1)
plot (etath1,I1,'-b','LineWidth',2)
hold on
plot(etath2,I2,'-r','LineWidth',2)
legend('Hydrogen @ M_\infty = 5','Methane @ M_\infty = 5')
legend('Location','north')
title('Specific impulse as a function of Thermal efficiency')
xlabel('\eta_t_h')
ylabel('I_s_p')



