function dydt=orbitfunperturbed(t,y)
%{
  y      - column vector containing the position and velocity vectors
           of the system at time t
  r      - position vector
  v      - velocity vector
  mu     - gravitational parameter for eartg
  rn     - magnitude of the relative position vector 
  a      - acceleration vectors of m1 & m2
  dydt   - column vector containing the velocity and acceleration
           vectors of the system at time t

%}
global mu pert
r=[y(1);y(2);y(3)];
v=[y(4);y(5);y(6)];
rn=norm(r);
a=[-(mu/rn)*r(1)+pert;-(mu/rn)*r(2)+pert;-(mu/rn)*r(3)+pert];
dydt=[v;a];
end