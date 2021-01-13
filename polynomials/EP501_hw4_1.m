%% EP 501 - Project 4
% Vittorio Baraldi
clear all; close all; clc

%% Exercise 1

%% Part a) and b)

load test_lsq.mat
addpath ../linear_algebra

[f,e,S]=leastsquares(x,ynoisy,1);
[f2,e2,S2]=leastsquares(x,ynoisy,2);
[f3,e3,S3]=leastsquares(x,ynoisy,3);

%checking results for cubic
poly=polyfit(x,ynoisy,3);

%outputs
figure(1)
plot(x,ynoisy,'-b','LineWidth',0.2)
hold on
plot(x,f,'--g','LineWidth',1.5)
hold on
plot(x,f2,'-r','LineWidth',1.5)
hold on
plot(x,f3,'-y','LineWidth',2)
hold on
plot(x,polyval(poly,x),'--k','LineWidth',1.2)
title('Linear Least Square Fit')
xlabel('x')
ylabel('y(x)')
legend('Noisy values','Linear fit','Quadratic Fit','Cubic Fit','MatLab polyfit for cubic')
figure(2)
plot(x,e,'-k','LineWidth',0.2)
hold on
plot(x,e2,'-m','LineWidth',0.2)
hold on
plot(x,e3,'-g','LineWidth',0.2)
title('Error comparation')
xlabel('x')
ylabel('Error')
legend('Linear fit error','Quadratic fit error','Cubic fit error')
disp('Linear fit residual:')
disp(S)
disp('Quadratic fit residual:')
disp(S2)
disp('Cubic fit residual:')
disp(S3)

%% Functions
%---------------------------------------------------------------------%

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

%------------------------------------------------------------------------%

function [Amod,ord]=Gauss_elim(A,b,verbose)

% [Amod,ord]=Gauss_elim(A,b,verbose)
%
% This function perform elimination with partial pivoting and scaling as
% described in Section 1.3.2 in the Hoffman textbook (viz. it does Gaussian
% elimination).  Note that the ordering which preserves upper triangularity
% is stored in the ord output variable, such that the upper triangular output
% is given by row-permuted matrix Amod(ord,:).  The verbose flag can be set to
% true or false (or omitted, default=false) in order to print out what the algirthm
% is doing for each elimination step.

%Parse the inputs, throw an error if something is obviously wrong with input data
narginchk(2,3);
if (nargin<3)
    verbose=false;
end %if

%Need to error check for square input.  

%Allocation of space and setup
Amod=cat(2,A,b);          %make a copy of A and modify with RHS of system
n=size(A,1);              %number of unknowns
ord=(1:n)';               %ord is a mapping from input row being operated upon to the actual row that represents in the matrix ordering

%Elimination with scaled, partial pivoting for matrix Amod; note all row
%indices must be screen through ord mapping.
for ir1=1:n-1
    if (verbose)
        disp('Starting Gauss elimination from row:  ');
        disp(ir1);
        disp('Current state of matrix:  ');
        disp(Amod(ord,:));
    end %if
    
    %check scaled pivot elements to see if reordering should be done
    pivmax=0;
    ipivmax=ir1;      %max pivot element should never be higher than my current position
    for ipiv=ir1:n    %look only below my current position in the matrix
        pivcurr=abs(Amod(ord(ipiv),ir1))/max(abs(Amod(ord(ipiv),:)));      %note that columns never get reordered...
        if (pivcurr>pivmax)
            pivmax=pivcurr;
            ipivmax=ipiv;     %this stores the index into ord for row having largest pivot element
        end %if
    end %for
    
    %reorder if situation calls for it
    if (ipivmax ~= ir1)
        itmp=ord(ir1);
        ord(ir1)=ord(ipivmax);
        ord(ipivmax)=itmp;
        
        if (verbose)
            disp('Interchanging rows:  ');
            disp(itmp);
            disp(' and:  ');
            disp(ord(ir1));
            disp('Current matrix state after interchange:  ');
            disp(Amod(ord,:));
        end %if
    end %if
    
    %perform the elimination for this row, former references to ir1 are now
    %mapped through the ord array
    for ir2=ir1+1:n
        fact=Amod(ord(ir2),ir1);
        Amod(ord(ir2),ir1:n+1)=Amod(ord(ir2),ir1:n+1)-fact/Amod(ord(ir1),ir1).*Amod(ord(ir1),ir1:n+1);    %only need columns ahead of where we are in matrix
    end %for
    
    if (verbose)
        disp('Following elimination for row:  ');
        disp(ir1);
        disp(' matrix state:  ');
        disp(Amod(ord,:));
    end %if
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