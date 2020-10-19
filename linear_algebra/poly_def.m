function g=poly_def(f,y)
syms x
%extracting coefficients from imported function
c=sym2poly(f(x));
c=flip(c);
n=length(c);
b=zeros(n,1);
b(n)=c(n);
for i=n-1:-1:1
    %calculating the new coefficient array b
    b(i)=c(i)+y*b(i+1);
end
%the first value of b is zero, therefore it does not apply
b=b(2:n);
b=flip(b);
%ply2sym is the opposite function of sym2poly, from coefficients given, it
%extracts the function in a symbolic function (sym 1x1)
g=poly2sym(b);
end
    