function [f,g,fdot,gdot] = LagrangeCoeffZ( r1m, r2m, dTheta, z, mu )

% Compute the Lagrange coefficients in terms of universal variable z
%
%   Inputs: 
%     r1m     Magnitude of position vector r1  
%     r2m     Magnitude of position vector r2
%     dTheta  Angle between vector r1 and r2
%     z       Universal variable "z"
%     mu      Gravitational constant 
% 
% 	Outputs:
%     f   
%     g
%     fdot
%     gdot 
%


C     = stumpC(z);
S     = stumpS(z);
A     = sin(dTheta)*sqrt(r1m*r2m/(1-cos(dTheta))); 
y     = r1m+r2m+A*( z*S-1 )/sqrt(C);                
f     = 1-y/r1m;                                    
g     = A*sqrt(y/mu);                              
fdot  = sqrt(mu)/(r1m*r2m)*sqrt(y/C)*(z*S-1);       
gdot  = 1-y/r2m;                                     