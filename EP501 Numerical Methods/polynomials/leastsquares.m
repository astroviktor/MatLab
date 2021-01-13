function [f,e,S]=leastsquares(x,y,n)
% This functions performs a least squares approximation of a polynomial of
% degree n
%{
Legend
n     -- degree of polynomial
x     -- x points
y     -- function noisy values
f     -- approximated function
e     -- error vector
S     -- final residual
%}
N=length(x);
%sumations of x^(ith) for matrix M
for i=1:2*n
    xi(i)=sum(x.^(i));
end
%rhs vector (b)
yi(1)=sum(y);
for i=2:n+1
    yi(i)=sum(x.^(i-1).*y);
end
%matrix M to solve for
M(1,1)=N;
for i=2:n+1
    M(1,i)=xi(i-1);
end
for i=2:n+1
    M(i,:)=[xi(i-1:n+i-1)];
end
%performing Gaussian Elimination and back-substitution
[Mmod,ord]=Gauss_elim(M,yi');
a=backsub(Mmod(ord,:));
%calculating the approximated function 
for i=1:N
    f(i)=a(1);
    for j=1:n
       f(i)=f(i)+a(j+1)*x(i)^j;
    end
        %error vector
       e(i)=f(i)-y(i);
end
%residual
S=e*e';

end