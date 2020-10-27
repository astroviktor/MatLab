function xv2=chisquared(f,y,sigma,n)

% This script evaluates the reduced chi-squared statistic, which defines a
% goodness-to-fit statistic of a polynomial of order n

%{
Legend
n = degree of polynomial
f = function values at points x_i
y = noisy values at points x_i
sigma = uncertainty of data
v = number of degrees of freedom in the fit
%}

v=length(y)-(n-1);
xv2=0;
for i=1:length(y)
    %evaluatin chi-sq. stat
    xv2=xv2+(((y(i)-f(i))^2)/sigma(i)^2);
end
xv2=xv2/v;
end