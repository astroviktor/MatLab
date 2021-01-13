%% EP 501 - Numerical Methods
%% Project 6
% Vittorio Baraldi
clear all; close all; clc;
addpath ../linear_algebra
%-------------------------------------------------------------------------%
%% Exercise 1

% Part a)
a=0.01;                 %[m]
l=a/5;                  %[m]
x1=-9*a/10;             %[m]
x2=9*a/10;              %[m]
e0=8.854e-12;           %[F/m]
x=linspace(-a,a,20);
n=length(x);
dx=x(2)-x(1);

%dielectric function
e=e0*(10*tanh((x-x1)/l)-10*tanh((x-x2)/l));

%dielectric function derivative
%forward difference at the beginning
de_dx(1)=(e(2)-e(1))/dx;
%centered difference on the interior
for i=2:n-1
    de_dx(i)=(e(i+1)-e(i-1))/2/dx;
end %for

%backward difference at the end
de_dx(n)=(e(n)-e(n-1))/dx;

%plots
figure;
plot(x,e,'LineWidth',2)
title('Dielectric function \epsilon(x)')
xlabel('x')
ylabel('\epsilon(x)')
figure;
plot(x,de_dx,'--k','LineWidth',2)
title('Dielectric function derivative')
xlabel('x')
ylabel('d\epsilon/dx')

% Part b) thru d)
% NOTE: System of equations for both electrostatic potential and boundary
% conditions in the next page.
A=zeros(n,n);
b=zeros(n,1);
%boundary conditions
b(1)=1000;
b(n)=100;

%using first order difference at
%1st point
A(1,1)=-1/dx; A(1,2)=1/dx;
A(n,n)=1;
for i=2:n-1
    for j=1:n       
        if i==j
            %diagonal points Y_i
            A(i,j)=-2*e(i)/dx/dx;
        elseif i==(j+1)
            %subdiagonal points Y_i-1
            A(i,j)=e(i)/dx/dx-de_dx(i)/2/dx;
        elseif i==(j-1)
            %superdiagonal points Y_i+1
            A(i,j)=e(i)/dx/dx+de_dx(i)/2/dx;
        end
    end
end

%solving system with built-in Matlab function (1st order difference)
sol_first=A\b;

% Part e)
%using second order difference
A2=A;
A2(1,1)=-3/2/dx; A2(1,2)=2/dx; A2(1,3)=-1/2/dx;

%solving system with built-in Matlab function (2nd order difference)
sol_second=A2\b;

%outputs and plots
figure;
plot(x,sol_first,'-b','LineWidth',1.2)
hold on
plot(x,sol_second,'--r','LineWidth',1.2)
title('Electrostatic Potential function \Phi(x)')
xlabel('x')
ylabel('\Phi(x)')
legend('First order difference','Second order difference')