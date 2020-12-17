%% EP 501 - Final
%% Vittorio Baraldi
%% EXERCISE 1
clear all; close all; clc;
addpath ../linear_algebra
addpath ../differentiation
addpath ../nonlinear_eqns
addpath ../polynomials

%% PART B)

%Amplification factor. The product alpha*delta-t is being considered as a
%one element in the function for simplicity. Alpha*delta-t is a positive
%value.
t=linspace(0,3,500);
tol=1e-10;
for i=1:length(t)
    %amplif. factor for RK4
    G(i)=1-t(i)+t(i)^2/2-t(i)^3/6+t(i)^4/24;
    %amplif. factor for RK2
    Grk2(i)=1-t(i)+t(i)^2/2;
    if t(i)>0
        if G(i)-1<tol
            flagx=t(i);
            flagy=G(i);
        end
        if Grk2(i)-1<tol
            flagrk2=t(i);
        end
    end
end

%outputs and plots
figure;
plot(t,G,'LineWidth',2)
hold on
xline(flagx,'--k')
hold on
scatter(flagx,flagy,70)
legend('Amplification Factor','Stability limit')
title('RK4 Method Stability')
xlabel('\alpha\Deltat')
ylabel('G(\alpha,\Deltat)')
set(gca,'FontSize',20)

%% PART C)

%Polynomials
funG1=@(x) -x+x^2/2-x^3/6+x^4/24; %G(alpha,dT)-1=0
funG2=@(x) 2-x+x^2/2-x^3/6+x^4/24;%G(alpha,dT)+1=0
xplot=linspace(0,4,200);
for j=1:length(xplot)
    g1(j)=funG1(xplot(j));
    g2(j)=funG2(xplot(j));
end
%maximum number of iterations
maxit=800;
%initial guess
x0=0.5+0.5i;
%all roots for function number 1 (funG1)
for i=1:4
    %finding a root with the approximated newton method (derivative
    %calculated numerically)
    [x1(i)]=newton_approx(funG1,x0,maxit,tol);
    %deflating the polynumial by the root found with the newton method
    funG1=poly_def(funG1,x1(i));
    %extracting the usable function from the sym 1x1 (i.e. you cannot assign
    %x,y... values to a sym 1x1)
    funG1=matlabFunction(funG1);
end

%all roots for function number 2 (funG2)
for i=1:4
    %finding a root with the approximated newton method (derivative
    %calculated numerically)
    [x2(i)]=newton_approx(funG2,x0,maxit,tol);
    %deflating the polynumial by the root found with the newton method
    funG2=poly_def(funG2,x2(i));
    %extracting the usable function from the sym 1x1 (i.e. you cannot assign
    %x,y... values to a sym 1x1)
    funG2=matlabFunction(funG2);
end
%%
%outputs and plots
figure;
plot(xplot,g1,'--r','LineWidth',1.5)
hold on
plot(xplot,g2,'--b','LineWidth',1.5)
title('|G(\alpha,\Deltat)| < 1')
xlabel('\alpha\Deltat')
ylabel('G(\alpha,\Deltat)')
legend('G(\alpha,\Deltat)-1=0','G(\alpha,\Deltat)+1=0')
set(gca,'FontSize',20)
figure;
plot(real(x1),imag(x1),'o','LineWidth',1.5,'MarkerSize',10)
hold on
plot(real(x2),imag(x2),'o','LineWidth',1.5,'MarkerSize',10)
hold on
xline(0)
hold on
yline(0)
legend('Roots for G-1=0','Roots for G+1=0')
title('Roots for |G(\alpha,\Deltat)| < 1')
xlabel('Real(\alpha\Deltat)')
ylabel('Imaginary(\alpha\Deltat)')
set(gca,'FontSize',20)
%%
% The only real solutions for the marginal stability limits are x=0 and
% x=2.7853


%% PART D)

% Find the greatest time step for RK4 and RK2 method
%amplif. factor for RK4
% Solving G-1=0, since it is the only option with real roots
Grk4=@(x) -x+x^2/2-x^3/6+x^4/24;
%amplif. factor for RK2
x0=2;
maxit=800;
tol=1e-10;
maxtimesteprk4=newton_approx(Grk4,x0,maxit,tol);
%%
% For the Runge Kutta Fourth-Order method, the greatest timestep that
% achieves stability is given by
%%
% $\Delta{t} = \frac
% {maxtimestep}{\alpha}$ 
%%
% while for the Second-Order method (RK2),
% the maximum timestep is given by
%%
% $\Delta{t}=
% \frac{1}{\alpha}$.
%%
% In the
% first expression, _maxtimestep_ is approximately equal to
disp(maxtimesteprk4)

%% PART E)
% Using a single-step evaluation for RK4 in this problem (see Part a)

% Picking a slightly larger value for the time step
adt=maxtimesteprk4*1.01;
t1=0:adt:30;
y0=1;
y1sol=RK4(t1,y0,@statescalar);
%Picking a slightly smaller value for the time step
adt=maxtimesteprk4*0.85;
t2=0:adt:30;
y2sol=RK4(t2,y0,@statescalar);

figure;
plot(t1,y1sol,'-r','LineWidth',1.5)
hold on
plot(t2,y2sol,'-b','LineWidth',1.5)
legend('Timestep greater than max','Timestep less than max')
title('ODE for different timesteps')
xlabel('Time')
ylabel('y(t)')
set(gca,'FontSize',15)

%% Functions 

function [root,it,success]=newton_approx(f,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f (input as a handle)
% given a function which computes the derivative
%Implements the fprime function numerically (derivative)

% Error checking of input
narginchk(2,5);   %check for correct number of inputs to function
if (nargin<3)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<4)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<5)
    verbose=false;
end %if
epsilon=.01;
% Newton iterations
it=1;
root=x0;
fval=f(root);
converged=false;
while(~converged && it<=maxit)
    fprime=(f(root+epsilon)-f(root))/epsilon; %derivative evaluation
    derivative=fprime;
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        root=root-fval./derivative;    % update root estimate
        fval=f(root);                  % see how far off we are from zero...
        if (verbose)
            fprintf(' iteration: %d; root:  %f + %f i; function value: %f, derivative:  %f \n',it,real(root),imag(root),fval,derivative);
        end %if
        it=it+1;
        converged=abs(fval)<tol;
    end %if
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function

%-------------------------------------------------------------------------%

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

%-------------------------------------------------------------------------%

function F=statescalar(t,U)
y=U(1);
F=-y;
end

%-------------------------------------------------------------------------%

function y=RK4(t,x0,f)
%This function solves an ODE using the Runge-Kutta Fourth-Order Method
n=length(t);
dt=t(2)-t(1);
[r,~]=size(x0);
y=[x0,zeros(r,n-1)];
for i=1:n-1
    % Partial updates
    dy1=f(t(i),y(:,i));
    dy2=f(t(i)+dt/2,y(:,i)+dy1*dt/2);
    dy3=f(t(i)+dt/2,y(:,i)+dy2*dt/2);
    dy4=f(t(i)+dt,y(:,i)+dy3*dt);
    % Function update
    y(:,i+1)=y(:,i)+dt*(dy1+2*dy2+2*dy3+dy4)/6;
end
end
