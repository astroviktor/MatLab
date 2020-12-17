%% Ex.3
clear all
close all

%% Part b
% initial conditions
theta_0=0;      
thetadot_0=0; 

tspan=0:0.05:10; 
g=9.81;
L=1;
w=sqrt(g/L);
a_d=cos(5*tspan)/pi/L;

options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y]=ode45(@bruteforceintegration,tspan,[theta_0;thetadot_0],options);

% plot brute force integral solution
figure(1)
plot(tspan,y(:,1),tspan,y(:,2),'LineWidth',2)
grid on
legend('\theta(t)','theta_dot(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Brute Force Integration Solution')

%% Part c

[t,c]=ode45(@variation,tspan,[theta_0;thetadot_0],options);

% use unperturbed solution with varying "initial conditions"
for i=1:size(c,1)
    theta_p(i)=c((i),1)*cos(w*tspan(i))+c((i),2)/w*sin(w*tspan(i)); % unperturbed position
    thetadot_p(i)=-w*c((i),1)*sin(w*tspan(i))+c((i),2)*cos(w*tspan(i)); % unperturbed veolcity
end

% plot unperturbed solution
figure(2)
plot(tspan,theta_p,tspan,thetadot_p,'LineWidth',2)
grid on
legend('\theta(t)','theta_dot(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Variation of Parameters Solution')

% Plot difference in methods
figure(3)
plot(tspan,theta_p'-y(:,1),tspan,thetadot_p'-y(:,2),'LineWidth',2)
grid on
legend('\theta(t)','theta_dot(t)')
xlabel('Time (sec)')
ylabel('Pos,Vel (m, m/s)')
title('Difference in Solutions')

%% Functions
function dx=bruteforceintegration(t,x)
g=9.81;
L=1;
w=sqrt(g/L);
a_d=cos(5*t)/pi/L;
dx(1)=x(2);
dx(2)=-w^2*x(1)-a_d;
dx=dx';
end

function dc=variation(t,c)
g=9.81;
L=1;
w=sqrt(g/L);
a_d=cos(5*t)/pi/L;
dc(1)=1/w*sin(w*t)*a_d;
dc(2)=-cos(w*t)*a_d;
dc=dc';
end