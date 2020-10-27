% Part f)

clear all; close all;clc

y=5/3;                          %adiabatic index
rho=1.67e-21;                   %plasma mass density [kg/m^3]
p=1.38e-11;                     %plasma thermal pressure [Pa]
B=10^(-9);                      %local magnetic field [T]
u0=1.256637062e-6;              %vacuum permeability [H/m]
Cs=sqrt(y*p/rho)/100000;        %adiabatic sound speed [m*10^-5/s]
Ca=sqrt(B^2/u0/rho)/100000;     %Alfvèn speed [m*10^-5/s]
theta=pi/4;                     %wave propagation angle [rad]
% wave function
f=@(v) (v^2-Ca^2*cos(theta)^2)*(Ca^2*Cs^2*cos(theta)^2-v^2*(Ca^2+Cs^2)+v^4);
% wave function derivative
fprime=@(v) 2*Ca^4*v*cos(theta)^2+4*Cs^2*Ca^2*v*cos(theta)^2-4*...
    Ca^2*v^3*cos(theta)^2-4*Ca^2*v^3-4*Cs^3*v^3+6*v^5;
%plotting the function to obtain initial guesses
v=linspace(-1.2,1.2,100);
for i=1:length(v)
    fg(i)=f(v(i));
end    
plot(v,fg)
grid on
xlabel('\nu')
ylabel('f(\nu)')
title('Wave function for \theta = \pi/4')
%initial guess
x0=Ca;
%max number of iterations
maxit=600;
%tolerance to achieve
tol=1e-9;
%repeated newton exact method for different initial points
[root1]=newton_exact(f,fprime,x0,maxit,tol);
x0=-Ca;
[root2]=newton_exact(f,fprime,x0,maxit,tol);
x0=13;
[root3]=newton_exact(f,fprime,x0,maxit,tol);
x0=-13;
[root4]=newton_exact(f,fprime,x0,maxit,tol);
x0=root1+0.0001;
[root5]=newton_exact(f,fprime,x0,maxit,tol);
x0=root2-0.0001;
[root6]=newton_exact(f,fprime,x0,maxit,tol);
sol=[root1;root2;root3;root4;root5;root6]*100000;
%outputs
for i=1:length(sol)
    fprintf('Root %i = %1.6f m/s\n',i,sol(i))
end
