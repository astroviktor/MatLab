%% EXERCISE 3
clear all; close all;clc
addpath ../linear_algebra
addpath ../differentiation
addpath ../nonlinear_eqns
addpath ../polynomials
addpath ../PDEs
%% PART B)

%Solving the system of equations derived in part a)
lambda=2;
dx=1/64;
dt=5*dx^2/2/lambda;
x=0:dx:1;
tmax=1024/((lambda*2*pi)/(2*dx)^2);
t=0:dt:tmax;
lt=numel(t);
lx=numel(x);
f=zeros(lx,lt);
f(:,1)=sin(2*pi*x)+sin(8*pi*x);

%allocate sparse array storage (this matrix is to be tridiag)
A=sparse(lx,lx);
b=zeros(lx,1);

%backward Euler-method in time with centered space differences
for n=1:lt-1
    A(1,1)=1;
    b(1)=0;
    for ix=2:lx-1
        %i-1 coeff
        A(ix,ix-1)=-lambda*dt/dx^2;
        
        %i coeff
        A(ix,ix)=1+2*lambda*dt/dx^2;
        
        %i+1 coeff
        A(ix,ix+1)=-lambda*dt/dx^2;
        
        b(ix)=f(ix,n);
    end %for
    A(lx,lx)=1;
    b(lx)=0;
    
    fnow=A\b;
    f(:,n+1)=fnow;
end %for

%outputs and plots
figure;
imagesc(t,x,f);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)')
title('First order BDF with CSD')
set(gca,'FontSize',15);

%% PART D)

%Solving the system of equations derived in part c)
lambda=2;
dx=1/64;
dt=5*dx^2/2/lambda;
x=0:dx:1;
tmax=1024/((lambda*2*pi)/(2*dx)^2);
t=0:dt:tmax;
lt=numel(t);
lx=numel(x);
f2=zeros(lx,lt);
f2(:,1)=sin(2*pi*x)+sin(8*pi*x);
f2(:,2)=f(:,2);

%allocate sparse array storage (this matrix is to be tridiag)
A=sparse(lx,lx);
b=zeros(lx,1);

%backward Euler-method in time with centered space differences
for n=3:lt
    A(1,1)=1;
    b(1)=0;
    for ix=2:lx-1
        %i-1 coeff
        A(ix,ix-1)=-lambda/dx^2;
        
        %i coeff
        A(ix,ix)=3/2/dt+2*lambda/dx^2;
        
        %i+1 coeff
        A(ix,ix+1)=-lambda/dx^2;
        
        b(ix)=2*f2(ix,n-1)/dt-f2(ix,n-2)/2/dt;
    end %for
    A(lx,lx)=1;
    b(lx)=0;
    
    fnow=A\b;
    f2(:,n)=fnow;
end %for
%plotting the second order BDF method
figure;
imagesc(t,x,f2);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)')
title('Second order BDF with CDS')
set(gca,'FontSize',15);

%exact solution
[T,X]=meshgrid(t,x);
exactsol=exp(-4*pi^2*lambda*T).*sin(2*pi*X)+exp(-64*pi^2*lambda*T)...
    .*sin(8*pi*X);

%Plotting exact solution
figure;
imagesc(t,x,exactsol);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('Exact solution for PDE');
set(gca,'FontSize',15);
        