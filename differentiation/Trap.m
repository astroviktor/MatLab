function [I] = Trap(fun,lmin,lmax,N)
%this function implements numerical integration via trapezoidal method


% For the simplest case, the calculation is done by the formula
if N==1
    I= 0.5*(lmax-lmin)*2*fun(1);
else

% For the multiple application, a sum is performed taking advantage of the 
% vector nature of the evaluation
    I=(0.5*(lmax-lmin)/N)*(fun(1)+fun(end)+2*sum(fun(2:end-1)));   

end
end







