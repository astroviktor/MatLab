%% TURBO RAMJET ANALYSIS
% Mixed cycle analysis (ByPass+Core)
clear all
close all
% Initial conditions
M0=[2];                       %mach number
g=1.4;                      %heat air coefficient
g0=9.81;                    %gravitational constant [m/s^2]
q=50000;                    %dynamic pressure [Pa]
R=287.058;                  %gas ideal constant [J/kg*K]
Ts=15;                      %atmospheric temperature at sea level [°C]
Ps=101325;                  %atmopsheric pressure at sea level [Pa]
q0=1+((g-1)/2)*M0.^2;
rhos=1.225;                 %atmospheric density at sea level [kg/m^3]
pt0_p0=q0.^(g/(g-1));

nd=.95;                     %diffuser efficiency
nc=.88;                     %compressor efficiency
nbp=.95;                    %bypass efficiency
nb=.96;                     %burner efficiency
nt=.92;                     %turbine efficiency
nab=.96;                    %after-burner efficiency
nn=.97;                     %nozzle efficiency
pi_c=12;                    %compressor pressure ratio
pt6_pt5=.94;
Tmax6=1550;                 %max. operative temperature for turbine [K]
pt10_pt9=.94;
Tt10=1900;                  %total temperature at the nozzle [K]
Q=45e6;                     %[J/kg]
Cp=1005;                    %heat capacity [J/kg*k]
Bc=.5;                      %fraction of mass flow through bypass
Bb=1-Bc;                    %fraction of mass flow through core
A0=1;                       %capture area [m^2]
A9=1;
A8=Bb;
A7=Bc;
Cf=0.005;                   %friction coefficient
L=1;                        %duct length
Aw=pi*(L/2)^2;          %area of the wall (1 meter lenght/1 meter diameter)
Ae=1;                       %exit area
for j=1:length(M0)
%% Altitude 

P0(j)=q*2/g/M0(j)^2;
%from the free-stream pressure we can consult tables and get the altitude
%and free-stream conditions
fprintf("Given the free-stream pressure P0=%i\n",P0(j))
rho0(j)=input("Insert density: ");
T0(j)=input("\nInsert temperature: ");
a0(j)=sqrt(g*R*T0(j));
u0(j)=M0(j)*a0(j);

Tt0(j)=T0(j)*q0(j);             %total free-stream temperature [k]
% Diffuser

pt3(j)=((1+nd*(q0(j)-1))^(g/(g-1)))*P0(j);
Tt3(j)=Tt0(j);

% Compressor

Tt5(j)=Tt3(j)*(1+(1/nc)*(pi_c^((g-1)/g)-1));
pt5(j)=pi_c*pt3(j);

% Burner

pt6(j)=pt6_pt5*pt5(j);
f(j)=(Tmax6-Tt5(j))/(nb*Q/Cp-Tmax6);

% Turbine

Tt7(j)=(Tt3(j)+Tmax6+f(j)*Tmax6-Tt5(j))/(1+f(j));
pt7(j)=((1-(1-Tt7(j)/Tmax6)/nt)^(g/(g-1)))*pt6(j);

%% Combined cycle
if Bc>0 && Bc<1
    %Mixer
    Tt8(j)=Tt0(j);
    pt8(j)=((1+nbp*(q0(j)-1))^(g/(g-1)))*P0(j);

    Tt9(j)=(Bc*Tt7(j)*(1+f(j))+Bb*Tt8(j))/((1+f(j))*Bc+Bb);
    % Here we guess a value for Mach number at point 7 and iterate. We use the
    % mass flow rate as check value. §Same applies for Mach at point 8
    % (following)
    M7g(j,1)=0.1;
    dm=0.0001;
    m0(j)=rho0(j)*u0(j)*A0;
    m7(j)=(1+f(j))*Bc*m0(j);
    for i=1:100000
        %calculating temperature, pressure and mass flow based on the guess
        q7g(j,i)=1+((g-1)/2)*M7g(j,i);
        T7g(j,i)=Tt7(j)/q7g(j,i);
        p7g(j,i)=pt7(j)/(q7g(j,i)^(g/(g-1)));
        m7g(j,i)=p7g(j,i)*sqrt(g*R*T7g(j,i))*M7g(j,i)*A7/(R*T7g(j,i));
        %check on guessed MFR and real one
        diff7(j,i)=m7(j)-m7g(j,i);
        l(j,i)=abs(diff7(j,i));
        if l(j,i)<10e-2
            p7(j)=p7g(j,i);
            T7(j)=T7g(j,i);
            M7(j)=M7g(j,i);
            break
        else
            M7g(j,i+1)=M7g(j,i)+dm;
        end
    end
    M8g(j,1)=.1;
    m8(j)=Bb*m0(j);
    for i=1:200000
        q8g(j,i)=1+((g-1)/2)*M8g(j,i);
        T8g(j,i)=Tt8(j)/q8g(j,i);
        p8g(j,i)=pt8(j)/(q8g(j,i)^(g/(g-1)));
        m8g(j,i)=p8g(j,i)*sqrt(g*R*T8g(j,i))*M8g(j,i)*A8/(R*T8g(j,i));
        diff8(j,i)=m8(j)-m8g(j,i);
        k(j,i)=abs(diff8(j,i));
        if k(j,i)<10e-2 %&& M8g(j,i)<1
            p8(j)=p8g(j,i);
            T8(j)=T8g(j,i);
            M8(j)=M8g(j,i);
            break
        else
            M8g(j,i+1)=M8g(j,i)+dm;
        end
    end
    % in order to get conditions at point 9, we need to get through a sort of
    % double iteration using the energy and momentum eqns. With them combined
    % we end up having an equation that we can solve for Mach 9. If such value
    % and our guess coincide we have found the correct answer.
    M9g(j,1)=.1;
    m9(j)=m7(j)+m8(j);
    syms P9
    for i=1:200000
        %momentum eqn to find pressure at point 9 using our guess
        num(j,i)=(p7(j)*A7*(g*M7(j)*M7(j)+1))+(p8(j)*A8*(g*M8(j)*M8(j)+1));
        den(j,i)=((g*M9g(j,i)*M9g(j,i)+1)*A9)+(Cf*g*Aw*M9g(j,i)^2)/2;
        p9g(j,i)=num(j,i)/den(j,i);
        q9g(j,i)=1+((g-1)/2)*M9g(j,i);
        T9g(j,i)=Tt9(j)/q9g(j,i);
        %using the energy eqn we find a new value for M9
        fac1(j,i)=((1+f(j))*Bc/Bb)+1;
        fac2(j,i)=(p8(j)/p9g(j,i))*(A8/A9)*(sqrt(T9g(j,i)/T8(j)));
        M9v(j,i)=fac1(j,i)*fac2(j,i)*M8(j);
        %checking on the guessed and calculated values to be equal (or close
        %enough!)
        diff9(j,i)=M9v(j,i)-M9g(j,i);
        y(j,i)=abs(diff9(j,i));
        if y(j,i)<10e-5 %&& M9g(j,i)<1
            M9(j)=M9g(j,i);
            p9(j)=p9g(j,i);
            T9(j)=T9g(j,i);
            q9(j)=q9g(j,i);
            break
        else
            M9g(j,i+1)=M9g(j,i)+dm;
        end
    end
    pt9(j)=p9(j)*(q9(j)^(g/(g-1)));

    %non-dimensional enthropy
    dS_Cp(j)=Bb*(log(T9(j)/T8(j))-((g-1)/g)*log(p9(j)/p8(j)))+(1+f(j))*Bc*...
        (log(T9(j)/T7(j))-((g-1)/g)*log(p9(j)/p7(j)));

    % After burner
    fab(j)=(1+f(j))*(Tt10-Tt9(j))/((nab*Q/Cp)-Tt10);

    % Nozzle
    pt10(j)=pt10_pt9*pt9(j);
    %exit velocity
    ue(j)=sqrt(2*nn*Cp*Tt10*(1-(P0(j)/pt10(j))^((g-1)/g)));
    %specific thrust
    T_m(j)=(1+f(j)+fab(j))*ue(j)-u0(j);
    %thrust specific fuel consumption
    TSFC(j)=(f(j)+fab(j))/T_m(j);

    %conditions at point 11 (exit)
    Tt11(j)=Tt10;
    p11_pt10(j)=P0(j)/pt10(j);
    T11(j)=Tt11(j)*(1-nn*(1-p11_pt10(j)^((g-1)/g)));
    a11(j)=sqrt(g*R*T11(j));
    M11(j)=ue(j)/a11(j);
    q11(j)=1+((g-1)/2)*M11(j)^2;
    P11(j)=P0(j);
    Pt11(j)=P11(j)*(q11(j)^(g/(g-1)));

   
%% Turbojet mode
elseif Bc==1
    %Mixer
    Tt9(j)=Tt7(j);
    pt9_pt0(j)=(pt7(j)/pt6(j))*(pt6(j)/pt5(j))*(pt5(j)/pt3(j))*...
        (pt3(j)/P0(j));
    pt10_pt0(j)=pt9_pt0(j)*pt10_pt9;
    %After burner
    fab(j)=(1+f(j))*(Tt10-Tt9(j))/((nab*Q/Cp)-Tt10);
    %Nozzle
    ue(j)=sqrt(2*nn*Cp*Tt10*(1-(1/pt10_pt0(j))^((g-1)/g)));
    T_m(j)=(1+f(j)+fab(j))*ue(j)-u0(j);
    TSFC(j)=(f(j)+fab(j))/T_m(j);
    Tt11(j)=Tt10(j);
    T11(j)=Tt11(j)*(1-nn*(1-(1/pt10_pt0(j))^((g-1)/g)));
    a11(j)=sqrt(g*R*T11(j));
    M11(j)=ue(j)/a11(j);
    q11(j)=1+((g-1)/2)*M11(j)^2;
    P11(j)=P0(j);
    Pt11(j)=P11(j)*(q11(j)^(g/(g-1)));
%% Ramjet mode
elseif Bc==0
    %Diffuser/Mixer
    Tt8(j)=Tt0(j);
    pt8(j)=(1+nbp*(q0(j)-1)^(g/(g-1)))*P0(j);
    Tt9(j)=Tt8(j);
    %After Burner (same pressure loss as for main burner in turbojet mode)
    fab(j)=(Tt10-Tt9(j))/((nab*Q/Cp)-Tt10);
    %Nozzle
   
    ue(j)=sqrt(2*nn*Cp*Tt10*(1-(P0(j)/pt8(j))^((g-1)/g)));
    T_m(j)=(1+fab(j))*ue(j)-u0(j);
    TSFC(j)=fab(j)/T_m(j);
    
    Tt11(j)=Tt10(j);
    p11_pt10(j)=P0(j)/pt8(j);
    T11(j)=Tt11(j)*(1-nn*(1-p11_pt10(j)^((g-1)/g)));
    a11(j)=sqrt(g*R*T11(j));
    M11(j)=ue(j)/a11(j);
    q11(j)=1+((g-1)/2)*M11(j)^2;
    P11(j)=P0(j);
    Pt11(j)=P11(j)*(q11(j)^(g/(g-1)));
end
end
%% Legend
%{
Tt# --> Total temperature
pt# --> total pressure
p# --> static pressure
T# --> static temperature
q --> dynamic pressure
q# --> Total/static conditions factor (based on mach no.)
m# --> mass flow rate
TSFC --> thrust specific fuel consumption
T_m --> specific thrust
u0/ue --> initial free-stream flight speed and exit velocity at the nozzle
          respectively
The sub_g means that we are considering a guessed value.
%}



