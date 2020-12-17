%% AE504 Adv. Compressible Flow
% Fanno Flow relations
function [k,ds,ratioP,ratioT,ratioR,ratioP0]=fannoflow(M,g)

%{
Legend
ratioP  -- static pressure/sonic pressure
ratioT  -- static temperature/sonic temperature
ratioR  -- static density/sonic density
ratioP0 -- total pressure/sonic tot. pressure
ds      -- entropy gain normalized to Cp (dS/Cp)
k       -- Fanno parameter (4*f*L/D)
M       -- Mach number
g       -- ratio of specific heat (gamma)

The function calculates the fanno parameter and the pressure, density and
temperature ratios for a given Mach number
%}
for i=1:length(M)
    if M(i)<=0
        fprintf('Mach no. at position %i is either null or negative',i)
    else
    k(i)=((1-M(i)^2)/g/M(i)^2)+((g+1)/2/g)*log(M(i)^2/...
        ((2/(g+1))*(1+((g-1)/2)*M(i)^2)));
    ds(i)=log((M(i)^((g-1)/g))*(((2/(g+1))*(1+((g-1)/2)*M(i)^2))^...
        ((-g-1)/2/g)));
    ratioP(i)=(1/M(i))*(1/sqrt((2/(g+1))*(1+((g-1)/2)*M(i)^2)));
    ratioR(i)=(1/M(i))*sqrt((2/(g+1))*(1+((g-1)/2)*M(i)^2));
    ratioT(i)=1/((2/(g+1))*(1+((g-1)/2)*M(i)^2));
    ratioP0(i)=(1/M(i))*(((2/(g+1))*(1+((g-1)/2)*M(i)^2))^...
        ((g+1)/(2*(g-1))));
    end

end
end