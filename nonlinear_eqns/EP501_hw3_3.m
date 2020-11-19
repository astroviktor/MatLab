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
%checking for very small values to approximate to 0, for output purposes
if rootx1<1e-5
    rootx1=0;
elseif rooty1<1e-5
    rooty1=0;
end
%second 2 roots
%initial guesses
x0=0;
y0=1.5;
[rootx2,rooty2]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol);
%checking for very small values to approximate to 0, for output purposes
if rootx2<1e-5
    rootx2=0;
elseif rooty2<1e-5
    rooty2=0;
end
x0=1.6-0.6i;
y0=-1-0.5i;
[rootx3,rooty3]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol);
x0=1.6+0.6i;
y0=-0.8+0.5i;
[rootx4,rooty4]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol);
%solution matrix: first column=x solutions, second column=y solutions
%outputs
fprintf("System roots:\n")
disp('x1:')
disp(rootx1)
disp('y1:')
disp(rooty1)
disp('x2:')
disp(rootx2)
disp('y2:')
disp(rooty2)
disp('x3:')
disp(rootx3)
disp('y3:')
disp(rooty3)
disp('x4:')
disp(rootx4)
disp('y4:')
disp(rooty4)
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
