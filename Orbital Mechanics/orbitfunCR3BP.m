function dydt=orbitfunCR3BP(~,y)
%{
  y      - column vector containing the position and velocity vectors
           of the system at time t
  r      - position vector
  v      - velocity vector
  mu     - gravitational parameter
  rn     - magnitude of the relative position vector 
  a      - acceleration vectors of m1 & m2
  dydt   - column vector containing the velocity and acceleration
           vectors of the system at time t

%}
global lambda 
r=[y(1);y(2)];
v=[y(3);y(4)];

%r1 and r2: norms
r1=sqrt((r(1)+lambda)^2+r(2)^2);
r2=sqrt((r(1)+lambda-1)^2+r(2)^2);

%expliciting acceleration terms in CR3BP EOMs
a1=((-(1-lambda)*(r(1)+lambda))/r1^3)-(lambda*(r(1)-1+lambda)/r2^3)+r(1)+2*v(2);
a2=((-(1-lambda)*r(2))/r1^3)-((lambda*r(2))/r2^3)+r(2)-2*v(1);

a=[a1;a2];
dydt=[v;a];
end