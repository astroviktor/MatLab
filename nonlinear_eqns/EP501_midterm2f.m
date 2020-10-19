% Part f)

clear all; close all;clc

y=5/3;
rho=1.67e-21;
p=1.38e-11;
B=10^(-9);
u0=1.256637062e-6;
Cs=sqrt(y*p/rho);
Ca=sqrt(B^2/u0/rho);
theta=pi/4;
f=@(v) -(Ca^2+Cs^2)*v^2+v^4+(Cs^2*Ca^2*cos(theta)^2);
fprime=@(v) -2*v*(Ca^2*Cs^2-2*v^2);
fprime2=@(v) -2*(Cs^2+Ca^2-6*v^2);
x1=Ca*cos(theta);
x2=-x1;
x0=8000;
maxit=7000;
tol=1e-10;
[root,it,success]=newton_approx(f,x0,maxit,tol);
[root,it,success]=newton_exact_multiple(f,fprime,fprime2,x0,maxit,tol);


