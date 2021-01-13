%% EP 501 Numerical Methods
%% Project 6
% Vittorio Baraldi
clear all;close all;clc;

%-------------------------------------------------------------------------%

%% Exercise 2
% Part a)
%Initial Data
m=1.67e-27;                     %[kg]
q=1.6e-19;                      %[C]
B=50000e-9;                     %[T]
T=2*pi*m/q/B;                   %oscillation period [s]
N=100;
t=linspace(0,T,N);
omega=q*B/m;                    %frequency of oscillation
dt=t(2)-t(1);                   %delta-t
lt=numel(t);
vx2=zeros(1,lt);
vy2=zeros(1,lt);
vx2(1)=1000;                     %vx initial condition [m/s]
vy2(1)=1000;                     %vy initial condition [m/s]

% RK 2 method
for n=2:lt
    %step x and y components together, this is the half update
    vxhalf=vx2(n-1)+dt/2*(omega*vy2(n-1));
    vyhalf=vy2(n-1)-dt/2*(omega*vx2(n-1));
    
    %now the full update
    vx2(n)=vx2(n-1)+dt*(omega*vyhalf);
    vy2(n)=vy2(n-1)-dt*(omega*vxhalf);    
end %for
%Integrate velocity to get position as a fn. of time, this assumes that the
%particle is initially at x,y = (0,0)
% RK 2
x2=cumtrapz(t,vx2);    %Matlab built-in for accumulating an integral value
y2=cumtrapz(t,vy2);
vz2=1000;               %z-velocity considered constant [m/s]
z2=vz2*t;                %z-coordinate of position is just velocity*time

% RK 4 Method          
initial=[0;0;vx2(1);vy2(1)];
sol=RK4(t,initial,@state);
vz=1000; 
z=vz*t;



%outputs and figures
figure;
comet3(x2,y2,z2)
set(gca,'FontSize',10);
xlabel('x');
ylabel('y');
zlabel('z');
title('Position as a function of time for RK2')

figure;
comet3(sol(1,:),sol(2,:),z)
set(gca,'FontSize',10);
xlabel('x');
ylabel('y');
zlabel('z');
title('Position as a function of time for RK4')

figure;
ax=plotyy(t,vx2,t,vy2);
set(ax(1),'FontSize',10);
set(ax(2),'FontSize',10);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Velocity components as functions of time for RK2')

figure;
ax=plotyy(t,sol(3,:),t,sol(4,:));
set(ax(1),'FontSize',10);
set(ax(2),'FontSize',10);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Velocity components as functions of time for RK4')

% Part b)
%Showing that RK4 can solve the problem with fewer time steps, i.e. it is
%more accurate than RK2

tb=linspace(0,T,2*N/3);            %half of the time steps than part a)
vxb=1000;                     %vx initial condition [m/s]
vyb=1000;                     %vy initial condition [m/s]
initialb=[0;0;vxb;vyb];
solb=RK4(tb,initialb,@state);
%integrals
% xb=cumtrapz(tb,vxb);
% yb=cumtrapz(tb,vyb);
vzb=1000; 
zb=vzb*tb;

%outputs and plots
figure;
ax=plotyy(tb,solb(3,:),tb,solb(4,:));
set(ax(1),'FontSize',10);
set(ax(2),'FontSize',10);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Velocity components as functions of time for RK4 (2/3 timesteps)')
figure;
comet3(solb(1,:),solb(2,:),zb)
set(gca,'FontSize',10);
xlabel('x');
ylabel('y');
zlabel('z');
title('Position as a function of time for RK4 (2/3 timesteps)')

% Part c)                   
tc=linspace(0,T*5,N*5);             %5 oscillations period
vxc=1000;                     %vx initial condition [m/s]
vyc=1000;                     %vy initial condition [m/s]
initialc=[0;0;vxc;vyc];
solc=RK4(tc,initialc,@state);

%figures and plots
figure;
ax=plotyy(tc,solc(3,:),tc,solc(4,:));
set(ax(1),'FontSize',10);
set(ax(2),'FontSize',10);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Velocity components as functions of time for RK4 (5 oscillations)')
figure;
plot(solc(1,:),solc(2,:))
set(gca,'FontSize',10);
xlabel('x(m)');
ylabel('y(m)');
title('Position as a function of time for RK4 (5 oscillations)')

%-------------------------------------------------------------------------%
%% Functions

% State Vector Function
function F=state(t,U)

q=1.6e-19;
m=1.67e-27;
B=50000e-9;
F=[U(3);U(4);q*(U(4)*B*(1+U(2)/2))/m;-q*(U(3)*B*(1+U(2)/2))/m];
end

%Runge-Kutta 4th order Method Function
function y=RK4(t,x0,f)
n=length(t);
dt=t(2)-t(1);
[j,~]=size(x0);
y=[x0,zeros(j,n-1)];
for i=1:n-1
    dy1=f(t(i),y(:,i));
    dy2=f(t(i)+dt/2,y(:,i)+dy1*dt/2);
    dy3=f(t(i)+dt/2,y(:,i)+dy2*dt/2);
    dy4=f(t(i)+dt,y(:,i)+dy3*dt);
    y(:,i+1)=y(:,i)+dt*(dy1+2*dy2+2*dy3+dy4)/6;
end
end