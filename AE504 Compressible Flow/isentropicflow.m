%% AE504 Adv. Compressible Flow
% Isentropic Flow Relations

function [ratioP,ratioT,ratioR,ratioTs,ratioPs,ratioRs]=isentropicflow(M,g);
%{ 
Legend
ratioP  -- static pressure/total pressure
ratioT  -- static temperature/total temperature
ratioR  -- static density/total density
ratioTs -- sonic temperature/total temperature
ratioPs -- sonic pressure/total pressure
ratioRs -- sonic density/total density
M       -- Mach number
g       -- heat air coefficent
%}
for i=1:length(M)
    if M(i)<=0
        fprintf('Mach no. at position %i is either null or negative',i)
    else   
        ratioP(i)=(1+((g-1)/2).*M(i)^2).^(-(g/(g-1)));
        ratioT(i)=ratioP(i)^((g-1)/g);
        ratioR(i)=ratioP(i)^(1/g);
        ratioTs(i)=2/(g+1);
        ratioPs(i)=ratioTs(i)^(g/(g-1));
        ratioRs(i)=ratioTs(i)^(1/(g-1));
    end
end
end
 