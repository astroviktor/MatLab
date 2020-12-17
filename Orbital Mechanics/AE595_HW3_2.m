%% ORBITAL MECHANICS 
% HOMEWORK 3

%% PART 2

a=12000; % semi major axis [km]
e=.25; %eccentricity
theta=120; %true anomaly [degrees]
mu=398600; %gravitational coefficient [km^3/s^2]
n=sqrt(mu/a^3); %mean orbit rate [rad/s]
thetarad=theta*pi/180;

t = TimeFromTrueAnom(e,thetarad,n); 

fprintf('Time of flight from perigee = %d seconds\n',t)

%% PART 3

a2=15000;
e2=.2;
t2=5100;
th0=0;
n2=sqrt(mu/a2^3);
h2=sqrt(a2*mu*(1-e2^2));

M2=n2*t2; %mean anomaly
Ef=@(x) x-e2*sin(x)-M2; %equation for Eccentric anomaly
dEf=@(x) 1-e2*cos(x);   %derivative of the Eccentric anomaly eqn

E2=newton_raphson(Ef,dEf,M2,1e-10,500)

truean=2*atan( sqrt((1+e2)./(1-e2)) .* tan(E2/2) ); %trueanomaly

r=(h2^2/mu)*(1/(1+e2*cos(truean))); %radius

fprintf('Radius after 5100 seconds = %d km\n',r)

