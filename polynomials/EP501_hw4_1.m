%% EP 501 - Project 4
% Vittorio Baraldi
clear all; close all; clc

%% Exercise 1

%% Part a) and b)

load test_lsq.mat
addpath ../linear_algebra
% sumations
xi=sum(x);
xi2=sum(x.^2);
xi3=sum(x.^3);
xi4=sum(x.^4);
xi5=sum(x.^5);
xi6=sum(x.^6);
yixi=sum(x.*ynoisy);
yixi2=sum(x.^2.*ynoisy);
yixi3=sum(x.^3.*ynoisy);
yi=sum(ynoisy);
%matrix form
M=[1*length(x) xi;xi xi2];
b=[yi;yixi];
%gaussian elimination + back substitution
[Mmod,ord]=Gauss_elim(M,b);
a=backsub(Mmod(ord,:));
%linear fit
for i=1:length(x)
    f(i)=a(1)+a(2)*x(i);
    e(i)=f(i)-ynoisy(i);
end
S=e*e';

M2=[1*length(x) xi xi2;xi xi2 xi3;xi2 xi3 xi4];
b2=[yi;yixi;yixi2];
[Mmod2,ord]=Gauss_elim(M2,b2);
a2=backsub(Mmod2(ord,:));
%quadratic fit
for i=1:length(x)
    f2(i)=a2(1)+a2(2)*x(i)+a2(3)*x(i)^2;
    e2(i)=f2(i)-ynoisy(i);
end
S2=e2*e2';

%cubic solution
M3=[1*length(x) xi xi2 xi3;xi xi2 xi3 xi4;xi2 xi3 xi4 xi5;xi3 xi4 xi5 xi6];
b3=[yi;yixi;yixi2;yixi3];
[Mmod3,ord]=Gauss_elim(M3,b3);
a3=backsub(Mmod3(ord,:));
%cubic fit
for i=1:length(x)
    f3(i)=a3(1)+a3(2)*x(i)+a3(3)*x(i)^2+a3(4)*x(i)^3;
    e3(i)=f3(i)-ynoisy(i);
end
S3=e3*e3';

%checking results for cubic
poly=polyfit(x,ynoisy,3);

%outputs
figure(1)
plot(x,ynoisy,'-b','LineWidth',0.2)
hold on
plot(x,f,'--g','LineWidth',1.5)
hold on
plot(x,f2,'-r','LineWidth',1.5)
hold on
plot(x,f3,'-y','LineWidth',2)
hold on
plot(x,polyval(poly,x),'--k','LineWidth',1.2)
title('Linear Least Square Fit')
xlabel('x')
ylabel('y(x)')
legend('Noisy values','Linear fit','Quadratic Fit','Cubic Fit','MatLab polyfit for quadratic')
figure(2)
plot(x,e,'-k','LineWidth',0.2)
hold on
plot(x,e2,'-m','LineWidth',0.2)
hold on
plot(x,e3,'-g','LineWidth',0.2)
title('Error comparation')
xlabel('x')
ylabel('Error')
legend('Linear fit error','Quadratic fit error','Cubic fit error')
disp('Linear fit residual:')
disp(S)
disp('Quadratic fit residual:')
disp(S2)
disp('Cubic fit residual:')
disp(S3)

%% Part c) and d)

stat=chisquared(f,ynoisy,sigmay,1);
stat2=chisquared(f2,ynoisy,sigmay,2);
stat3=chisquared(f3,ynoisy,sigmay,3);



disp('According to the statistic chi-squared approach, for the above function we can see how the cubic method works best (i.e. closest value to 1)')
disp('Chi-squared values:')
disp('n=1')
disp(stat)
disp('n=2')
disp(stat2)
disp('n=3')
disp(stat3)