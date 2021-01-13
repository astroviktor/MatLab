function x=tridiag(A,b)

%This function solves a tridiagonal system of equations using the Thomas
%algorithm
% A = matrix (needs to be tridiagonal)
% b = rhs vector

n=length(b);
k=b;
A1(1,1)=0;
A1(n,3)=0;
%creating the A' matrix
    %first column
    for i=2:n
        A1(i,1)=A(i,i-1);
    end
    %second column
    for i=1:n
        A1(i,2)=A(i,i);
    end
    %third column
    for i=1:n-1
        A1(i,3)=A(i,i+1);
    end
       
    %matrix A' and vector b elimination   
    for i=2:n
        A1(i,2)=A1(i,2)-(A1(i,1)/A1(i-1,2))*A1(i-1,3);
        k(i)=k(i)-(A1(i,1)/A1(i-1,2))*k(i-1);
    end
    
    x=zeros(n,1);
    %back substitution
    for i=n-1:-1:1
        x(n)=k(n)/A1(n,2);
        x(i)=(k(i)-A1(i,3)*x(i+1))/A1(i,2);
    end
end