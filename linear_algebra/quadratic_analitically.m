function [x1,x2]=quadratic_analitically(f)
syms x
%extracting coefficients from imported function
c=sym2poly(f(x));
a=c(1);
b=c(2);
c=c(3);
%calculating delta
delta=b^2-4*a*c;
    if delta<0
        % impossible equation
        disp('The equation does not admit solutions...')
    elseif delta==0
        % undetermined eqn
        x1=-b/2/a;
        x2=x1;
    else
        %two real roots
        x1=(-b+sqrt(delta))/2/a;
        x2=(-b-sqrt(delta))/2/a;
    end
end