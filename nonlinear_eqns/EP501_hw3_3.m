%% EP 501 - Numerical Methods
%% Homework 3
% Vittorio Baraldi

clear all; close all; clc

%% Ex. 3
%% Part a)

%functions
f=@objfun2Df;
g=@objfun2Dg;
%functions' gradients
gradf=@grad_objfun2Df;
gradg=@grad_objfun2Dg;
%first 2 roots
%initial guesses
x0=1.5;
y0=0;
maxit=1000;  %max number of iterations allowed
tol=1e-10;   %tolerance to achieve
[rootx1,rooty1,it]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol);
%second 2 roots
%initial guesses
x0=0;
y0=1.5;
[rootx2,rooty2]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol);
%solution matrix: first column=x solutions, second column=y solutions
sol=[rootx1 rooty1;rootx2 rooty2];

%cleaning values smaller than 10^-8 (for output purposes)
for i=1:size(sol,1)
    for j=1:size(sol,1)
    if abs(sol(i,j))<1e-8
        sol(i,j)=0;
    end
    end
end
%outputs
fprintf("System roots:\n")
fprintf('x=%1.2f\t',sol(1,1))
fprintf('y=%1.2f\n\n',sol(1,2))
fprintf('x=%1.2f\t',sol(2,1))
fprintf('y=%1.2f\n\n',sol(2,2))

%% Part b)

%functions (3D problem)
f2=@objfun3Df;
g2=@objfun3Dg;
h=@objfun3Dh;
%functions' gradients
gradf2=@grad_objfun3Df;
gradg2=@grad_objfun3Dg;
gradh=@grad_objfun3Dh;
%initial guesses
x0=0.5;
y0=1.2;
z0=2;
tol2=1e-10;
maxit2=600;
[rootx,rooty,rootz]=newton3D_exact(f2,gradf2,g2,gradg2,h,gradh,x0,y0,z0,maxit2,tol2);

%outputs
fprintf('3-D System root;\n')
fprintf('x = %1.4f\ty = %1.4f\tz = %1.4f\n',rootx,rooty,rootz)
