function dydt=orbitfun(t,y,mu)
%{
  y      - column vector containing the position and velocity vectors
           of the system at time t
  r      - position vector
  v      - velocity vector
  mu     - gravitational parameter for earth
  rn     - magnitude of the relative position vector 
  a      - acceleration vectors of m1 & m2
  dydt   - column vector containing the velocity and acceleration
           vectors of the system at time t

%}
r=[y(1);y(2);y(3)];
v=[y(4);y(5);y(6)];

rn=norm(r);
%Ceres' mass is sufficiently big, but in this case we do not want to assume
%the spacecraft mass to be neglectable for the equation of motion.
G=6.67430e-20;
%Dawn's mission launch mass
m_space=1217.7; 
m_ceres=9.3835e20;
avec=G*(m_space+m_ceres)/(rn^3);
%we can consider the orbit to be non-perturbed since in the asteroid belt,
%both Jupiter and Mars are far away enough and do not give a solid
%contribution. Moreover, other bodies in the belt are either very small
%and/or distant from Ceres so for a 6 months orbit they do not apply.
a=[avec*r(1);avec*r(2);avec*r(3)];

dydt=[v;a];
end