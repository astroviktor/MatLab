%% EP 501 - Numerical Methods
%% Project 5
clc; clear all; close all
% Vittorio Baraldi
%% Exercise 1
%Part a)
%initial data
mu0=4*pi*10^(-7);             %[H/m]
I=10;                      %[A]
a=.005;                    %[m]
x=linspace(-3*a,3*a,1000);
y=linspace(-3*a,3*a,1000);[xx,yy]=meshgrid(x,y);
for j=1:length(x)
    for i=1:length(x)
    %condition checking
    if sqrt(xx(i,j)^2+yy(i,j)^2)<a
        %x and y component of magnetic field
        Bx(i,j)= (mu0*I/2/pi/a/a)*sqrt(xx(i,j)^2+yy(i,j)^2)*(-yy(i,j)/sqrt(xx(i,j)^2+yy(i,j)^2));
        By(i,j)= (mu0*I/2/pi/a/a)*sqrt(xx(i,j)^2+yy(i,j)^2)*(xx(i,j)/sqrt(xx(i,j)^2+yy(i,j)^2));
    elseif sqrt(xx(i,j)^2+yy(i,j)^2)>=a
        Bx(i,j)= (mu0*I/2/pi)/sqrt(xx(i,j)^2+yy(i,j)^2)*(-yy(i,j)/sqrt(xx(i,j)^2+yy(i,j)^2));
        By(i,j)= (mu0*I/2/pi)/sqrt(xx(i,j)^2+yy(i,j)^2)*(xx(i,j)/sqrt(xx(i,j)^2+yy(i,j)^2));
    end
    end
end
%outputs
figure(1)
pcolor(xx,yy,Bx)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{x}')
colorbar
figure(2)
pcolor(xx,yy,By)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{y}')
colorbar

%part b)
x2=linspace(-3*a,3*a,20);
y2=linspace(-3*a,3*a,20);[xx2,yy2]=meshgrid(x2,y2);
for j=1:length(x2)
for i=1:length(x2)
    %condition checking
    if sqrt(xx2(i,j)^2+yy2(i,j)^2)<a
        %x and y component of magnetic field
        Bx2(i,j)= (mu0*I/2/pi/a)*sqrt(xx2(i,j)^2+yy2(i,j)^2)*(-yy2(i,j)/sqrt(xx2(i,j)^2+yy2(i,j)^2));
        By2(i,j)= (mu0*I/2/pi/a)*sqrt(xx2(i,j)^2+yy2(i,j)^2)*(xx2(i,j)/sqrt(xx2(i,j)^2+yy2(i,j)^2));
    elseif sqrt(xx2(i)^2+yy2(i)^2)>=a
        Bx2(i,j)= (mu0*I/2/pi)/sqrt(xx2(i,j)^2+yy2(i,j)^2)*(-yy2(i,j)/sqrt(xx2(i,j)^2+yy2(i,j)^2));
        By2(i,j)= (mu0*I/2/pi)/sqrt(xx2(i,j)^2+yy2(i,j)^2)*(xx2(i,j)/sqrt(xx2(i,j)^2+yy2(i,j)^2));
    end
end
end
%outputs
figure(3)
quiver(xx2,yy2,Bx2,By2)
title('Quiver plot')
xlabel('x')
ylabel('y')

%Part c)

dx=x(2)-x(1);
dy=y(2)-y(1);
gradx=zeros(size(Bx));
grady=zeros(size(By));

%x component of gradient of By
for i=2:length(x)-1
    gradx(:,i)=(By(:,i+1)-By(:,i-1))/2/dx; 
end
gradx(:,1)=(By(:,2)-By(:,1))/dx;
gradx(:,length(x))=(By(:,length(x))-By(:,length(x)-1))/dx;
%y component of gradient Bx
for i=2:length(y)-1
    grady(:,i)=(Bx(:,i+1)-Bx(:,i-1))/2/dy; 
end
grady(:,1)=(Bx(:,2)-Bx(:,1))/dy;
grady(:,length(y))=(Bx(:,length(y))-Bx(:,length(y)-1))/dy;
%calculating curl (z-component only)
curlz=gradx-grady;
%outputs
figure(4)
imagesc(x,y,curlz)
shading flat
xlabel('x')
ylabel('y')
title('Plot for curl of vector field B (numerical)')
colorbar

%part d)
%Derivatives computed with Wolfram-Alpha
for j=1:length(x)
for i=1:length(x)
    %condition checking
    if sqrt(xx(i,j)^2+yy(i,j)^2)<a
        %x and y component of magnetic field
        Bx3(i,j)= -0.08;
        By3(i,j)= 0.08;
    elseif sqrt(xx(i,j)^2+yy(i,j)^2)>=a
        Bx3(i,j)= (yy(i,j)^2-xx(i,j)^2)/500000/(xx(i,j)^2+yy(i,j)^2)^2;
        By3(i,j)= Bx3(i,j);
    end
end
end
curlzhand=-Bx3+By3;
figure(5)
imagesc(x,y,curlzhand)
shading flat
xlabel('x')
ylabel('y')
title('Plot for curl of vector field B (hand-computed)')
colorbar

%% Ex. 4
% Part a)

phi=linspace(0,2*pi,1000); dphi=phi(2)-phi(1); lx=length(phi);
phimesh=meshgrid(phi);
r0=2*a;
rx=r0.*cos(phimesh);
ry=r0.*sin(phimesh);
r=rx+ry;
%outputs
figure(6)
pcolor(xx,yy,Bx)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{x}')
colorbar

figure(7)
pcolor(xx,yy,By)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{y}')
colorbar

% Part b)
for j=1:length(x)
    for i=1:length(x)
    %condition checking
    if sqrt(rx(i,j)^2+ry(i,j)^2)<a
        %x and y component of magnetic field
        Bx4(i,j)= (mu0*I/2/pi/a/a)*sqrt(rx(i,j)^2+ry(i,j)^2)*(-ry(i,j)/sqrt(rx(i,j)^2+ry(i,j)^2));
        By4(i,j)= (mu0*I/2/pi/a/a)*sqrt(rx(i,j)^2+ry(i,j)^2)*(rx(i,j)/sqrt(rx(i,j)^2+ry(i,j)^2));
    elseif sqrt(rx(i,j)^2+ry(i,j)^2)>=a
        Bx4(i,j)= (mu0*I/2/pi)/sqrt(rx(i,j)^2+ry(i,j)^2)*(-ry(i,j)/sqrt(rx(i,j)^2+ry(i,j)^2));
        By4(i,j)= (mu0*I/2/pi)/sqrt(rx(i,j)^2+ry(i,j)^2)*(rx(i,j)/sqrt(rx(i,j)^2+ry(i,j)^2));
    end
    end
end
%outputs
figure(8)
pcolor(xx,yy,Bx4)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{x}')
colorbar

figure(9)
pcolor(xx,yy,By4)
shading flat
xlabel('x')
ylabel('y')
title('Plot for B_{y}')
colorbar

% Part c)

%Numerical Derivative
%centered difference on the interior
for j=1:lx
    %forward difference at the beginning
    dx_dphi(j,1)=(rx(j,2)-rx(j,1))/dphi;
    dy_dphi(j,1)=(ry(j,2)-ry(j,1))/dphi;   
    for ix=2:lx-1
        dx_dphi(j,ix)=(rx(j,ix+1)-rx(j,ix-1))/2/dphi;
        dy_dphi(j,ix)=(ry(j,ix+1)-ry(j,ix-1))/2/dphi;
    end %for
    %backward difference at the end
    dx_dphi(j,lx)=(rx(j,lx)-rx(j,lx-1))/dphi;
    dy_dphi(j,lx)=(ry(j,lx)-ry(j,lx-1))/dphi;
end
%Hand-computed derivative
drx=-r0.*sin(phimesh);
dry=r0.*cos(phimesh);
%summing components
dr=drx+dry;
dr_dphi=dx_dphi+dy_dphi;
%outputs and plots
figure(10)
plot(phi,dr_dphi(1,1:end),'-b','LineWidth',2)
hold on
plot(phi,dr(1,1:end),'--r','LineWidth',1.5)
title('Derivative of r(\phi) path')
xlabel('\phi')
ylabel('\phi')
zlabel('r(\phi)')
legend('Numerical derivative','Hand-computed derivative')

% Part d)

B=Bx4+By4;
f=B*dr_dphi/mu0;
I=0;
for i=1:length(phi)
    I=I+Trap(f,0,2*pi,1000);
end
disp('Auxiliary Magnetic Field (integral value):')
disp(I)