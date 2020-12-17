%% Phobos' Gravitational Potential
% Author: Vittorio Baraldi
% Contact: baraldiv@my.erau.edu

%------------------------------------------------------------------------%
clear all; close all; clc;
global mu R_p M pert
%------------------------------------------------------------------------%

[l,m,Cnorm,Snorm,C,S]=readvars('Phobos_gravitational_coefficients.txt');
for j=1:length(l)
    Cnorm2(j-m(j),m(j)+1)=Cnorm(j);
    Snorm2(j-m(j),m(j)+1)=Snorm(j); 
end
Cnorm=Cnorm2(any(Cnorm2,2),:);
Snorm=[zeros(1,11);Snorm2(any(Snorm2,2),:)];
%reference radius [km] (from reference paper)
R_p=14;
%reference mass [kg] (from reference paper)
M=1.06e16;
%gravitational parameter for Phobos [km^3/s^2]
mu=7.113588120963051e-4;
e=0.01;
a=100;
W=pi/20;
w=0;
i=1*pi/180;
vec=[Cnorm;Snorm];
th=0;
[r,v] = RVFromCOE( a,i,W,w,e,th, mu );
rn=norm(r);
vn=norm(v);
longitude=asin(r(3)/rn);
latitude=atan2(r(2),r(1));

for i=1:10
    P=legendre(i,sin(latitude));
end
for i=1:11
    for j=1:11
        U(i,j)=(((R_p/(R_p+rn))^(i-1))*P(i)*...
         (Cnorm(i,j)*cos(m(j)*longitude)+Snorm(i,j)*sin(m(j)*longitude)));
    end
end
U=U*(mu/(R_p+rn));
pert=-gradient(U);
pert=sum(sum(pert));
    y0=[r;v];
    tspan=[0 2000];
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,yf]=ode113(@orbitfunperturbed,tspan,y0,options);
    figure(1)
    plot3(yf(:,1),yf(:,2),yf(:,3))

