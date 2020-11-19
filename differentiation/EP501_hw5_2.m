%% EP 501 - Numerical Methods
%% Project 5
% Vittorio Baraldi
clear all;close all;clc;
%% Ex. 2
% Part a)
Q=1;            %[C]
a=1;            %[m]
e0=8.854e-12;   %[F/m]
x=linspace(-3*a,3*a,800);
y=linspace(-3*a,3*a,800);[xx,yy]=meshgrid(x,y);
%calculating the scalar field
for j=1:length(x)
    for i=1:length(x)
        %checking conditions
        check(i,j)=sqrt(xx(i,j)^2+yy(i,j)^2);
        if check(i,j)<a
            phi(i,j)=Q/4/pi/e0/a-Q/8/pi/a^3/e0*(xx(i,j)^2+yy(i,j)^2-a^2);
        elseif check(i,j)>=a
            phi(i,j)=Q/4/pi/e0/check(i,j);
        end
    end
end

%outputs
figure(1)
pcolor(xx,yy,phi)
shading flat
xlabel('x')
ylabel('y')
title('Plot for \Phi(x,y)')
colorbar

% Part b)
% Computing Laplacian
% First we need to compute the gradient
gradx=zeros(size(phi));
grady=zeros(size(phi));
dx=x(2)-x(1);
dy=y(2)-y(1);
%taking the gradient
%x component of gradient
for i=2:length(x)-1
    gradx(:,i)=(phi(:,i+1)-phi(:,i-1))/2/dx;    
end 
gradx(:,1)=(phi(:,2)-phi(:,1))/dx;
gradx(:,length(x))=(phi(:,length(x))-phi(:,length(x)-1))/dx;

%y component of gradient
for i=2:length(y)-1
    grady(i,:)=(phi(i+1,:)-phi(i-1,:))/2/dy;    
end 
grady(1,:)=(phi(2,:)-phi(1,:))/dy;
grady(length(y),:)=(phi(length(y),:)-phi(length(y)-1,:))/dy;

%taking the Laplacian
f=gradx;
g=grady;

%x derivative part of the divergence
divx=zeros(size(f));
for i=2:length(x)-1
    divx(:,i)=(f(:,i+1)-f(:,i-1))/2/dx;
end %for
divx(:,1)=(f(:,2)-f(:,1))/dx;
divx(:,length(x))=(f(:,length(x))-f(:,length(x)-1))/dx;

%y derivative part of the divergence
divy=zeros(size(y));
for i=2:length(y)-1
    divy(i,:)=(g(i+1,:)-g(i-1,:))/2/dy;
end %for
divy(1,:)=(g(2,:)-g(1,:))/dy;
divy(length(y),:)=(g(length(y),:)-g(length(y)-1,:))/dy;

%summing components
div=divx+divy;

%outputs and figures
figure(2)
surface(x,y,div)
grid on
shading flat
colorbar
xlabel('x')
ylabel('y')
title('Plot for numerical Laplacian of \Phi(x,y)')

% Part c)

% computing gradient by hand
% Derivatives calculated via Wolfram-Alpha
for j=1:length(x)
    for i=1:length(x)
    check(i,j)=sqrt(xx(i,j)^2+yy(i,j)^2);
    if check(i,j)<a
        divphix(i,j)=-8996887682;
        divphiy(i,j)=-8996887682;
    elseif check(i,j)>=a
        divphix(i,j)=8987742438*(3*xx(i,j)^2/(xx(i,j)^2+yy(i,j)^2)^(5/2))-(1/(xx(i,j)^2+yy(i,j)^2)^(3/2));
        divphiy(i,j)=8987742438*(3*yy(i,j)^2/(xx(i,j)^2+yy(i,j)^2)^(5/2))-(1/(xx(i,j)^2+yy(i,j)^2)^(3/2));

    end
    end
end

%summing components
divphi=divphix+divphiy;

%outputs and plots
figure(3)
surface(x,y,divphi)
grid on
shading flat
colorbar
xlabel('x')
ylabel('y')
title('Plot for hand-computed Laplacian of \Phi(x,y)')

%% Ex.3
z=linspace(-3*a,3*a,800);
zz=meshgrid(z);

%calculating phi and laplacian with z-coordinate included
for i=1:length(x) 
    for j=1:length(x)
        check(i,j)=sqrt(xx(i,j)^2+yy(i,j)^2+zz(i,j)^2);
        if check(i,j)<a
            div3(i)=-8996887682;
            phi3(i)=Q/4/pi/e0/a-Q/8/pi/a^3/e0*(xx(i,j)^2+yy(i,j)^2+zz(i,j)^2-a^2);
        elseif check(i,j)>=a
            phi3(i,j)=Q/4/pi/e0/check(i,j);
            num1=3*yy(i,j)^2;
            exp=(xx(i,j)^2+yy(i,j)^2+zz(i,j)^2)^(5/2);
            rr=num1/exp;
            %Divergence is equal in expression for each partial derivation with
            %respect the variables, plus all the variables have the same
            %length, same boundary values and equally spaced
            div3(i,j)=8987742438*rr-(1/(xx(i,j)^2+yy(i,j)^2+zz(i,j)^2)^(3/2));
        end
        f(i,j)=(e0*div3(i,j))*phi(i,j);
    end   
end
I=0;
for i=1:length(x)
    %computing integral via Trap subfunction (triple integral)
    I=I+Trap(Trap(Trap(f(i,:),-3*a,3*a,800),-3*a,3*a,800),-3*a,3*a,800);
end
%computing W_e (electrostatic energy)
We=-0.5*I;
disp('Electrostatic energy value (We):')
disp(We)

