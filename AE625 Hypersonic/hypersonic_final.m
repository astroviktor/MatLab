%% AE 625 - Final Exam
% Vittorio Baraldi

close all
clear all
%% Ex. 1
T=1000;                 %temperature [k]
Tr=2.9;                 %rotationmal temperature [K]
k=1.38064852e-23;       %Boltzmann's constant [J/K]
h=6.62607015e-34;       %Planck's constant [J*s]
beta=1/k/T;
I=h^2/(8*k*Tr*pi^2);    %moment of inertia for a diatomic nitrogen molecule
                        %(from literature) [kg*m^2]
                        
Q=T/Tr;                 %partition function derived in class
E=k*T;                  %internal rotational energy derived in class

%since index cannot have a zero/non-integer value, we calculate epsilon(0)
%and Q_rot(0) separately
g_0=1;
l_0=0;
epsilon_0=Tr*k*(l_0+1)*l_0;
Q_rot(1)=(2*l_0+1)*exp(-l_0*(l_0+1)*Tr/T);
Nj_N0=(1/Q)*g_0*exp(-epsilon_0/k/T);
for l=1:10000
    %energy level at j
    epsilon(l)=l*(l+1)*k*Tr;
    %degeneracy at j
    g(l)=2*l+1;
    %partition function calculated using sumation 
    Q_rot(l+1)=Q_rot(l)+((2*l+1)*exp(-l*(l+1)*Tr/T));
    %N_j/N
    Nj_N(l)=(1/Q)*g(l)*exp(-epsilon(l)/k/T);
    %checking on values. If both epsilon and the sum of Q_rot's are equal
    %or over the 99% of the values derived in class, we break the loop
    ratioQ(l)=Q_rot(l+1)/Q;
    ratioE(l)=epsilon(l)/E;
    if ratioQ(l)>=.99 && ratioE(l)>=.99
        x=l;
        break
    end
end

%% Ex. 2
th_r=2.5;               %rotational characteristic temperature for NO [K]
th_vib=2740;            %vibrational characteristic temperature for NO [K]

T2=[1000:1000:6000];    %temperature array
R=8.314;                %gas constant [J/mol*K]

for i=1:length(T2)
    %translation energy
    e_trans(i)=3*R*T2(i)/2;
    %vibration energy
    e_vib(i)=R*th_vib/(exp(th_vib/T2(i))-1);
    %rotation energy
    e_rot(i)=R*T2(i);
    
    %Cv's are given by the partial derivative of the corresponding energy
    %with respect to temperature at constant volume
    
    %Rotation Cv
    cvrot(i)=R;
    %Translation Cv
    cvtrans(i)=3*R/2;
    %Vibration Cv
    cvvib(i)=(th_vib^2*R*exp(th_vib/T2(i)))/(T2(i)^2*...
        (exp(th_vib/T2(i))-1)^2);
    %Summing the Cv's contribution from every energy
    cv(i)=cvvib(i)+cvrot(i)+cvtrans(i);
    %Given that R=Cp-Cv, we can obtain Cp 
    cp(i)=cv(i)+R;
    %The heat capacity ratio is given by the total Cp over the total Cv
    gamma(i)=cp(i)/cv(i);
end

%% Ex. 3
th_vib3=2270;           %vibrational characteristic temperature for O2 [K]
T3=[1000:1000:5000];    %temperature array [K]

for i=1:length(T3)
    %Calculating the entropy of vibration per molecule (derivation is in
    %the document
    S_N(i)=k*(-log(1-exp(-th_vib3/T3(i)))+(th_vib3/T3(i))/...
        (exp(th_vib3/T3(i))-1));
end