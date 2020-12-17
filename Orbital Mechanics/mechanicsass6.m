%Mechanics of Flight
%Assignment 6
%Problem 1.2

t0=518.67; %temperature in °R at sea level
a0=-0.003567; %rate of temperature decreasing in troposphere
h=[20000:-10:0]; %energy height vector in ft
R=1716.49; %gas constant
gamma=1.4; 
g=32.2;
v=linspace(0,800,length(h)); %velocity vector generation (ft/s)

for i=1:length(h)
   h0(i)=h(i)-(v(i)^2)/(2*g); %actual altitude
   t(i)=t0+a0*h0(i);%temperature at a given altitude
   c(i)=sqrt(gamma*R*t(i)); %speed of sound at a given altitude
   M(i)=v(i)/c(i); %mach number
   
end

plot(M,h0,'-b')
hold on
plot(M(501),h0(501),'r*')
grid on
xlabel('Mach no.')
ylabel('Altitude (ft)')
