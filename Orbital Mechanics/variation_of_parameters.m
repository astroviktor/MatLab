% Example of variation of parameters
% using spring-mass with perturbation

clear all; close all; clc;

% initial conditions
x_0=1;       % position, m
xdot_0=0.25; % velocity, m/s

t=0:0.05:10; % time span of sim, seconds

w=1;         % frequency, w^2=k/m

%% unperturbed solution
% x_u=x_0*cos(w*t)+xdot_0/w*sin(w*t); % unperturbed position
% xdot_u=-w*x_0*sin(w*t)+xdot_0*cos(w*t); % unperturbed veolcity
% 
% % plot unperturbed solution
% figure
% plot(t,x_u,t,xdot_u,'LineWidth',2)
% grid on
% legend('x_u(t)','xdot_u(t)')
% xlabel('Time (sec)')
% ylabel('Pos,Vel (m, m/s)')
% title('Unperturbed Solution')

%% brute force integration of perturbed problem
tspan=t;
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[T,X]=ode45(@brute,tspan,[x_0;xdot_0],options);

% plot brute force integral solution
figure
plot(t,X(:,1),t,X(:,2),'LineWidth',2)
grid on
legend('x_p(t)','xdot_p(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Brute Force Integration Solution')

%% variation of parameters solution
tspan=t;
% options=odeset('RelTol',1e-9,'AbsTol',1e-9);
[T,C]=ode45(@vop,tspan,[x_0;xdot_0],options);

% use unperturbed solution with varying "initial conditions"
for i=1:size(C,1)
    x_p(i)=C((i),1)*cos(w*t(i))+C((i),2)/w*sin(w*t(i)); % unperturbed position
    xdot_p(i)=-w*C((i),1)*sin(w*t(i))+C((i),2)*cos(w*t(i)); % unperturbed veolcity
end

% plot unperturbed solution
figure
plot(t,x_p,t,xdot_p,'LineWidth',2)
grid on
legend('x_p(t)','xdot_p(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Variation of Parameters Solution')

%% plot difference in methods
figure
plot(t,x_p'-X(:,1),t,xdot_p'-X(:,2),'LineWidth',2)
grid on
legend('x_p(t)','xdot_p(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Difference in Solutions')
% Notice that the difference is the value of options, line 27

%% FUNCTIONS
function dx=brute(t,x)

% let the perturbation be a_d=exp(-0.1*w*t)
w=1;

dx(1)=x(2);
dx(2)=-w^2*x(1)+exp(-.1*w*t);

dx=dx';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dc=vop(t,c)

% let the perturbation be a_d=exp(-0.1*w*t)
w=1;
a_d=exp(-0.1*w*t);

dc(1)=-1/w*sin(w*t)*a_d;
dc(2)=cos(w*t)*a_d;

dc=dc';
end