function [x,it]=newton_raphson(f,j,x0,toll,maxit)
% Use the Newton-Rhapson method to find the zero crossing for a function.
% 
%   Inputs:
%     f     Function handle. We want to find x where fun(x) = 0.
%     j     Derivative function handle. 
%     x0    Initial guess.
%     toll  Tolerance
%     maxit Maximum number of iteration (to avoid loops)
%   Outputs:
%     x         Solution to fun(x) = 0.
%     it        number of iterations
arrest=1;
it=0;

while arrest>toll && norm(f(x0))>=toll && it<=maxit
    s0=-j(x0)\f(x0);
    x=x0+s0;  
    arrest=norm(x-x0)/norm(x);
    x0=x;
    it=it+1;
end
it=it-1;
