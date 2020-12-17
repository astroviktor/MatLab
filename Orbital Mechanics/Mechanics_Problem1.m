%Problem 1

S=47;
AR=6.5;
e=0.87;
W=103047;
Ta=40298*2; %thrust available at sea level
rho=0.7364; %density at 5km of altitude
cd0=0.032;
rhosl=1.225; %density at sea level
%generating a velocity vector

V=[100:30:310];

for i=1:length(V)
    %calculating drag and lift coefficient
    cl(i)=(2*W)/(rho*V(i)^2*S);
    k=1/(pi*e*AR);
    cdi(i)=(k^2)*cl(i);
    cd(i)=cd0+cdi(i);
    %calculating thrust required and power required
    Tr(i)=W/(cl(i)/cd(i));
    Pr(i)=Tr(i)*V(i);
end

%calculatin power available from engines at sea level
%provided the thrust eqn in pont b) of the problem
Ta=Ta*(rho/rhosl);

for i=1:length(V)
    Pa(i)=Ta*V(i);
end

%plotting results
plot(V,Pr,'r');
grid on
hold on
plot(V,Pa,'b');
xlabel('Velocity');
ylabel('Power');
title('Power Required/Available at 5km');
legend('Power required', 'Power available');


