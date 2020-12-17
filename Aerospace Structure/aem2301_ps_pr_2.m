clear all;
close all;
clc;
%  Time and density history for the launch and ascent of a rocket
t = [0 10 20 30 40 50 60 70 80 90 100];                     %
rho = (1e-3)*[2.377 2.208 2.010 2.0624 2.1149 2.1673 2.2197 2.2721 2.3246 2.3770 2.377];
%   Loop through t and density to calculate altitude
number_of_data_points = length(t);
for k=1:number_of_data_points
  %   Density in slugs/ft^3
       h(k) = rho2h(rho(k));  
   end
   %   Plot figures
   figure(1)
   plot(t,rho,'r-');grid on;xlabel('Time');ylabel('Density (slugs/ft^3)');
   figure(2)
   plot(t,h,'r-');grid on;xlabel('Time');ylabel('Altitude (ft)');