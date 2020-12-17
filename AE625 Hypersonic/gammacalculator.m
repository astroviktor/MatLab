function gamma=gammacalculator(T)
n=length(T);
Cv=zeros(1,n);
gamma=zeros(1,n);
R=8.3145 %gas constant (J/mol*K)

%here follow the functions for the single energy components of a gas (x
%represents the temperature)

%etrans= 3/2*R*x; -- internal transition energy
%erot= R*x; -- internal rotational energy
%evibox=(R*2270)/(exp(2270/x)-1); -- internal vibration energy for oxygen
%evibn=(R*3390)/(exp(3390/x)-1); -- internal vibration energy for nitrogen
%eelec=(R*11390)*(2/3*exp(-11390/x))/(1+(2/3*exp(-11390/x))); -- electronic energy

%derivatives of the above functions - the purpose of derivate them is to
%directly obatin the Cv when summing up the computed energy states
%devibox=(5152900*exp(2270/x)*R)/(((-1+exp(2270/x))^2)*x^2);
%devibn=(11492100*exp(3390/x)*R)/(((-1+exp(3390/x))^2)*x^2);
%deelec=(778392600*exp(11390/x)*R)/(((2+(3*exp(11390/x)))^2)*x^2);

detrans=3/2*R;
derot=R;        %Cv for both the transition and rotational energy states is constant
devibox=zeros(1,n);
devibn=zeros(1,n);
deelec=zeros(1,n);

for i=1:n
    devibox(i)=0;
    devibn(i)=0;
    deelec(i)=0;
    if T(i)>=2270 %condition for oxygen
        devibox(i)=(5152900*exp(2270/T(i))*R)/(((-1+exp(2270/T(i)))^2)*T(i)^2);
        if T(i)>=3390 %condition for nitrogen
            devibn(i)=(11492100*exp(3390/T(i))*R)/(((-1+exp(3390/T(i)))^2)*T(i)^2);
        end
    end
    if T(i)>=11390 %condition for oxygen; nitrogen does not affect
        deelec(i)=(778392600*exp(11390/T(i))*R)/(((2+(3*exp(11390/T(i))))^2)*T(i)^2); %electronic energy
    end
    Cv(i)=detrans+derot+devibox(i)+devibn(i)+deelec(i); %total Cv
    gamma(i)=1+(R/Cv(i)); %heat capacity ratio
end
%--plots--
figure(1)
plot (T,gamma)
xlabel('Temperature')
ylabel('Heat Capacity ratio')
title('HC ratio with respect to temperature (for air)')

%contributes of the electronic and vibration energy states
figure(2)
plot (T,deelec,T,devibox,'--',T,devibn,'-r')
xlabel('Temperature')
ylabel('Energy contribution')
title('Contribute per energy state')
    
    
