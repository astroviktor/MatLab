function [value, isterminal, direction] = odecheck(~,y)
global lambda a b
% value = (abs(y(2)-land(2))+abs(y(1)-land(1)) < tol);

x = y(1); y = y(2);
% renaming 

value = (x-1+lambda)^2/a^2 + y^2/b^2 - 1; 
% checking if s/c is "inside ellipse," i.e. inside Phobos

isterminal = 1;   % Stop the integration
direction  = 0;
end