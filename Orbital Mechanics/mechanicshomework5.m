%Mechanics of Flight
%Homework 5
%Problem 4

TW=0.3797; %thrust to weight ratio
WS=76.84; %weight to span ratio
cd0=0.015; %parasite drag coeff
K=0.08; %induced drag coeff factor
clmax=1.2; %max lift coefficient

rho=2.377e-3; %density at sea level

v=[100:10:1250]; %generating velocity vector

nmaxcl=0.5.*rho.*v.^2.*(clmax/WS); %max load factor (constraint about clmax)
nmaxt=sqrt(((0.5.*rho.*v.^2)/(K.*WS)).*(TW-0.5.*rho.*v.^2.*(cd0/WS))); %max load factor (constraint about thrust)

nmaxtot=zeros(1,length(v));
for i=1:length(v)
    nmaxtot(i)=max(nmaxcl(i),nmaxt(i)); %calculating the actual max load factor based on the comparison between the two constraints
end

%plot
plot(v,nmaxt,':r',v,nmaxcl,':g',v,nmaxtot,'-.b')
grid on
xlabel('Velocity (ft/s)')
ylabel('Load Factor')
title('Load Factor vs Velocity')
legend('Load Factor (thrust)','Load Factor (stall)','Load Factor (actual)')




