function M=n_order_taylor(n)
% This function generates a matrix of coefficients to solve the system of
% equation to obtain an arbitrary derivative of order n

%checking for even numbers
if mod(n,2)==0
    for k=-n/2:n/2
        i=k+1+n/2;
        for j=1:n
            % taylor series expansion coefficients for each derivative
            M(i,j)=((k)^j)/factorial(j);
        end
    end
else
    n=n+1;
    for k=-n/2:(n/2-1)
        i=k+1+n/2;
        for j=1:n-1
            M(i,j)=((k)^j)/factorial(j);
        end
    end
end
% eliminating the k=0 terms
M=M(any(M,2),:);
% we need the inverse matrix in order to solve for the derivative
M=inv(M);
end