function g=poly_def(f,y)
%this function performes polynomial deflation
%{
Legend:
g=deflated polynomial
f=initial polynomial
y=polynomial's root
%}

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
%the first value of b is zero, so it does not apply
b=b(2:n);
b=flip(b);
%converting coefficient into function
g=poly2sym(b);
end