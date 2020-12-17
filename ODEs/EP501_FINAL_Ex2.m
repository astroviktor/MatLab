%% EXERCISE 2
clear all; close all;clc
addpath ../linear_algebra
addpath ../differentiation
addpath ../nonlinear_eqns
addpath ../polynomials

%% PART B)
% Testing the third-derivative expression obtained via Taylor-Series
% expansion:
%%
% $f'''(x_i)=\frac{3f(x_i)-f(x_{i-1})+
% f(x_{i+2})-3f(x_{i+1})}{\Delta{x}^3}$

x=linspace(0,2*pi,100);
dx=x(2)-x(1);
f=cos(x);
%evaluating derivative at center points
for i=2:length(x)-2
    f3(i)=(3*f(i)-3*f(i+1)+f(i+2)-f(i-1))/dx^3;
end
%outputs and plots
figure;
plot(x(2:end-1),f3,'-r')
hold on
plot(x,sin(x),'--b')
legend('Numerical solution','Exact solution')
title('Third-derivative of f=cos(x)')
xlabel('x')
ylabel('sin(x)')
set(gca,'FontSize',15)

%% PART E)
% Numerically invert the matrix M in the system
%%
% $f' = M^{-1}\Delta{f}$

%Generating the inverse of matrix M via LU Factorization
M=[-2 +2 -4/3 +2/3; -1 1/2 -1/6 1/24;...
    1 1/2 1/6 1/24; 2 2 4/3 2/3];
[L,U,bt]=simple_elimination_dolittle(M,[1;0;0;0]);
n=size(M,1);
for i=1:n
    %initializing the RHS
    bt=zeros(n,1);
    %setting the i-th element to 1, all others to 0
    bt(i)=1;
    %performing a forwards+back substitutions
    inv1=forwsub([L bt]);
    inv2=backsub([U inv1]);
    %inv2 is equal to the i-th column of the inverse matrix
    invmatrix(:,i)=inv2;
end
%checking result with buil-in MatLab function
checkmatrix=inv(M);
disp('Inverse matrix via LU Factorization:')
disp(invmatrix)
disp('Inverse matrix via MatLab built-in function:')
disp(checkmatrix)

%% PART F)
% Computing the fourth derivative derived in part _f_
for i=3:length(x)-2
    f4(i)=(-4*(f(i-1)+f(i+1))+f(i-2)+f(i+2)+6*f(i))/dx^4;
end
% The fourth derivative of f=cos(x) is equal to f itself

%outputs and plots
figure;
plot(x(3:end-2),f4(3:end),'-r','LineWidth',1.5)
hold on
plot(x,f,'--b','LineWidth',1.5)
legend('Numerical','Hand-computed')
xlabel('x')
ylabel('f^{(4)}(x)')
title('Fourth-derivative of f=cos(x)')
set(gca,'FontSize',15)

%% PART G)
% Writing a script that delivers the $M^{-1}$ matrix for a derivative or
% arbitrary order n

M5=n_order_taylor(5);
M6=n_order_taylor(6);
%%
disp('Matrix for 5th derivative')
disp(M5)
disp('Matrix for 6th derivative')
disp(M6)

%%
% The formulas for the 5th and 6th derivative of f=cos(x) are
%%
% $f^{(5)}_i=
% \frac{-f_{i-3}+5f_{i-2}-10f_{i-1}-
% 5f_{i+1}+f_{i+2}+10f_i}{\Delta{x}^5}$
%%
% $f^{(6)}_i=
% \frac{f_{i-3}-6f_{i-2}+15f_{i-1}+15f_{i+1}
% -6f_{i+2}+f_{i+3}-20f_i}{\Delta{x}^6}$

%% Functions

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

%-------------------------------------------------------------------------%

function [L,U,x]=simple_elimination_dolittle(A,b)
%This function provides a simple forward elimination method as already
%implemented in class examples that can be used with any matrix A and any
%vector, using the dolittle method

nref=length(b);    %system size for reference problem
L=eye(nref);
%note that the elimination procedure coded below modifies the matrix B
Awork=cat(2,A,b);   %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
for ir1=2:nref       %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref    %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=Awork(ir2,ir1-1)/Awork(ir1-1,ir1-1);     %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        Awork(ir2,:)=Awork(ir2,:)-fact.*Awork(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements), this is a little bit wasteful as it uses entire row...
        L(ir2,ir1-1)=fact; %lower triangular matrix
    
    end %for
end %for
x=Awork(:,nref+1);      %final solution
U=Awork(1:nref,1:nref); %upper triangular matrix
end

%-------------------------------------------------------------------------%

function x=forwsub(A)

% This function performs forward substitution on a lower triangular matrix


n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(1)=A(1,n+1)/A(1,1);          %finalized solution for last variable, resulting from upper triangular conversion

for ir1=2:n
    x(ir1)=A(ir1,n+1);       %assume we're only dealing with a single right-hand side here.
    fact=A(ir1,ir1);         %diagonal element to be divided through doing subs for the ir2 row
    for ic=ir1-1:-1:1
        x(ir1)=x(ir1)-A(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for

end %function

%-------------------------------------------------------------------------%

function x=backsub(A)

% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.


n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);          %finalized solution for last variable, resulting from upper triangular conversion

for ir1=n-1:-1:1
    x(ir1)=A(ir1,n+1);       %assume we're only dealing with a single right-hand side here.
    fact=A(ir1,ir1);         %diagonal element to be divided through doing subs for the ir2 row
    for ic=ir1+1:n
        x(ir1)=x(ir1)-A(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for

end %function
