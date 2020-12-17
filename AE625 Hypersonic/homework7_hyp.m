%% AE 625 - HOMEWORK 7
%% EX. 4.1
clear all
close all

t = [1:0.1:20];   %array of initial T_3/T_0
g = 1.36;         %heat air coefficient
syms m01 m02
for i=1:length(t)
    %solving eq. 4.9 pag. 158 for M3=0 and M3=1
    M0_1(i)=sqrt((t(i)-1)*(2/(g-1)));
    
    M0_2(i)=sqrt(((t(i)*(((g-1)/2)+1)-1)*(2/(g-1))));
    
end

plot(t,M0_1,'-r','LineWidth',1.2)
hold on
plot(t,M0_2,'-b','LineWidth',1.2)
ylabel('M_0')
xlabel('{T_3}/_{T_0}')
legend('M_3 = 0','M_3 = 1','location','north')
title('Free stream Mach number as a function of ^{T_3}/_{T_0}')
    
