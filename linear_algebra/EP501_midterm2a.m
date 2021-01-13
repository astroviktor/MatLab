%% Exercise 3

close all; clear all; clc
syms x
% Part a)
%quadratic function
fa=@(x) 2*x^2-6*x+4;
%matlab function for quadratic function resolution
[xa1,xa2]=quadratic_analitically(fa);
%outputs
disp('Solutions:')
disp(xa1)
disp(xa2)

% Part b)
%initial function
fb=@(x) x^5-15*x^4+85*x^3-225*x^2+274*x-120;
%1st solution (given)
x1=5;
%plynomial deflation
g=poly_def(fb,x1);
%outputs
disp('Plynomial divided by (x-5):')
disp(g)

% Part c) thru d)
%initial guess
x0=0.5;
%maximum number of iterations allowed for newton method
maxit=500;
%tolerance to achieve
tol=1e-10;
i=1;
%check parameter
converged=false;
%check function (equal to 1)
check=@()1.0;
while ~converged 
    %finding a root with the approximated newton method (derivative
    %calculated numerically)
    [x(i)]=newton_approx(fb,x0,maxit,tol);
    fprintf('Root number %i:\t%1.4f\n',i,x(i))
    %deflating the polynumial by the root found with the newton method
    fb=poly_def(fb,x(i));
    %extracting the usable function from the sym 1x1 (i.e. you cannot assign
    %x,y... values to a sym 1x1)
    fb=matlabFunction(fb);
    %once all the roots are found, the deflated function should be equal to
    %1.0
    converged=strcmp(func2str(check),func2str(fb));
    i=i+1;
end