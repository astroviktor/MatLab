%% ORBITAL MECHANICS - MIDTERM
%% Ex. 8
clear all
close all
%initial conditions
a1=12000;           %semi major axis [km]
e1=.2;              %eccentricity
i1=45*pi/180;       %inclination [rad]
w1=0;               %argument of perigee [rad]
omega1=0;           %RAAN [rad]
f1=0;               %true anomaly [rad]
p1=a1*(1-e1^2);

%final conditions
a2=13000;
e2=.1;
i2=45*pi/180;
w2=0;
omega2=0;
f2=120*pi/180;
p2=a2*(1-e2^2);

dt=60*60;           %transfer time [sec]
mu=398600;          %gravitational parameter [km^3/s^2]

%position and velocity vectors from orbital elements
[r1vec,v1vec] = RVFromCOE( a1,i1,omega1,w1,e1,f1, mu );
[r2vec,v2vec] = RVFromCOE( a2,i2,omega2,w2,e2,f2, mu );

%calculation of chord vector
cvec=r2vec-r1vec;

r2=norm(r2vec);
r1=norm(r1vec);

%angle between the two position vectors
theta=acos(r1vec'*r2vec/r1/r2);
c=sqrt(r1^2+r2^2-2*r1*r2*cos(theta));
%parameter "s"
s=(c+r1+r2)/2;
%picking max and min values for semi major axis of the transfer orbit
amin=s/2;
amax=100000*amin;

%minimum energy ellipse values
betamin=2*asin(sqrt((s-c)/s));
tmin=sqrt((amin^3/mu)*(pi-betamin+sin(betamin)));

%iterations to obtain the actual semi major axis of the transfer orbit
t=tmin;
g=t-dt;
while abs(g)>1e-10
    a=(amin+amax)/2;
    alpha=2*asin(sqrt(s/(2*a)));
    beta=2*asin(sqrt((r1+r2-c)/(4*a)));
    t=sqrt((a^3/mu)*(alpha-sin(alpha)-(beta-sin(beta))));
    if t<dt
        amin=a;
    elseif t>dt
        amax=a;
    end
    g=t-dt;
end

%solving for deltaE and p given the semi major axis found with the
%iteration (according to the procedure in the "Fundamentals of
%Astrodynamics" notes)
syms x y
eqn1=1-(r2/x)*(1-cos(theta))==1-(a/r1)*(1-cos(y));
eqn2=(r1*r2*sin(theta))/sqrt(mu*x)==t-sqrt(a^3/mu)*(y-sin(y));

sol = vpasolve([eqn1, eqn2], [x, y]);
xSol = double(sol.x);
ySol = double(sol.y);

%lagrange coefficients
f2=1-(a/r1)*(1-cos(ySol));
g2=t-sqrt(a^3/mu)*(ySol-sin(ySol));
fdot2=sqrt(mu/xSol)*tan(theta/2)*((1-cos(theta))/xSol-1/r1-1/r2);
gdot2=1-(r1/xSol)*(1-cos(theta));

%velocity at point 1 and 2 for the transfer orbit
vf1=1/g2*(r2vec-f2*r1vec);
vf2=1/g2*(gdot2*r2vec-r1vec);

%delta v's
dv1=norm(vf1-v1vec);
dv2=norm(v2vec-vf2);

%Orbital elements for the transfer orbit
[at,it,omegat,wt,et] = OrbitalElementsFromRV( r1vec, vf1, mu );

fprintf('deltaV 1 = %.2f km/s\ndeltaV 2 = %.2f km/s\n\n',dv1,dv2)
fprintf('Orbital elements of the transfer orbit\n')
fprintf('Semi major axis: %.0f km\n',at)
fprintf('Inclination: %.0f°\n',it*180/pi)
fprintf('Right ascension of the ascending node: %.0f°\n',omegat*180/pi)
fprintf('Argument of periapsis: %.1f°\n',wt*180/pi)
fprintf('Eccentricity: %.2f\n',et)    

    
   