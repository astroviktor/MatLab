clc;
close all;
clear all;
gamma=1.4;
M=1.01:0.01:7 
for i=1:600
    theeta_1=5;
    theeta_2=10;
    theeta_3=15;
    theeta_4=20;
    theeta_5=30;
    theta_1=theeta_1*pi/180;
    theta_2=theeta_2*pi/180;
    theta_3=theeta_3*pi/180;
    theta_4=theeta_4*pi/180;
    theta_5=theeta_5*pi/180;% convert to radians
    mu(i)=asin(1/M(i));                   % Mach wave angle
    c(i)=tan(mu(i))^2;
    a_1(i)=((gamma-1)/2+(gamma+1)*c(i)/2)*tan(theta_1);
    b_1(i)=((gamma+1)/2+(gamma+3)*c(i)/2)*tan(theta_1);
    a_2(i)=((gamma-1)/2+(gamma+1)*c(i)/2)*tan(theta_2);
    b_2(i)=((gamma+1)/2+(gamma+3)*c(i)/2)*tan(theta_2);
    a_3(i)=((gamma-1)/2+(gamma+1)*c(i)/2)*tan(theta_3);
    b_3(i)=((gamma+1)/2+(gamma+3)*c(i)/2)*tan(theta_3);
    a_4(i)=((gamma-1)/2+(gamma+1)*c(i)/2)*tan(theta_4);
    b_4(i)=((gamma+1)/2+(gamma+3)*c(i)/2)*tan(theta_4);
    a_5(i)=((gamma-1)/2+(gamma+1)*c(i)/2)*tan(theta_5);
    b_5(i)=((gamma+1)/2+(gamma+3)*c(i)/2)*tan(theta_5);
    d_1(i)=sqrt(4*(1-3*a_1(i)*b_1(i))^3/((27*a_1(i)^2*c(i)+9*a_1(i)*b_1(i)-2)^2)-1);
    beeta_1(i)=atan((b_1(i)+9*a_1(i)*c(i))/(2*(1-3*a_1(i)*b_1(i)))-(d_1(i)*(27*a_1(i)^2*c(i)+9*a_1(i)*b_1(i)-2))/(6*a_1(i)*(1-3*a_1(i)*b_1(i)))*tan(0*pi/3+1/3*atan(1/d_1(i))))*180/pi;
    d_2(i)=sqrt(4*(1-3*a_2(i)*b_2(i))^3/((27*a_2(i)^2*c(i)+9*a_2(i)*b_2(i)-2)^2)-1);
    beeta_2(i)=atan((b_2(i)+9*a_2(i)*c(i))/(2*(1-3*a_2(i)*b_2(i)))-(d_2(i)*(27*a_2(i)^2*c(i)+9*a_2(i)*b_2(i)-2))/(6*a_2(i)*(1-3*a_2(i)*b_2(i)))*tan(0*pi/3+1/3*atan(1/d_2(i))))*180/pi;
    d_3(i)=sqrt(4*(1-3*a_3(i)*b_3(i))^3/((27*a_3(i)^2*c(i)+9*a_3(i)*b_3(i)-2)^2)-1);
    beeta_3(i)=atan((b_3(i)+9*a_3(i)*c(i))/(2*(1-3*a_3(i)*b_3(i)))-(d_3(i)*(27*a_3(i)^2*c(i)+9*a_3(i)*b_3(i)-2))/(6*a_3(i)*(1-3*a_3(i)*b_3(i)))*tan(0*pi/3+1/3*atan(1/d_3(i))))*180/pi;
    d_4(i)=sqrt(4*(1-3*a_4(i)*b_4(i))^3/((27*a_4(i)^2*c(i)+9*a_4(i)*b_4(i)-2)^2)-1);
    beeta_4(i)=atan((b_4(i)+9*a_4(i)*c(i))/(2*(1-3*a_4(i)*b_4(i)))-(d_4(i)*(27*a_4(i)^2*c(i)+9*a_4(i)*b_4(i)-2))/(6*a_4(i)*(1-3*a_4(i)*b_4(i)))*tan(0*pi/3+1/3*atan(1/d_4(i))))*180/pi;
    d_5(i)=sqrt(4*(1-3*a_5(i)*b_5(i))^3/((27*a_5(i)^2*c(i)+9*a_5(i)*b_5(i)-2)^2)-1);
    beeta_5(i)=atan((b_5(i)+9*a_5(i)*c(i))/(2*(1-3*a_5(i)*b_5(i)))-(d_5(i)*(27*a_5(i)^2*c(i)+9*a_5(i)*b_5(i)-2))/(6*a_5(i)*(1-3*a_5(i)*b_5(i)))*tan(0*pi/3+1/3*atan(1/d_5(i))))*180/pi;
    M_n_1(i)=M(i)*sind(beeta_1(i));
    M_n_2(i)=M(i)*sind(beeta_2(i));
    M_n_3(i)=M(i)*sind(beeta_3(i));
    M_n_4(i)=M(i)*sind(beeta_4(i));
    M_n_5(i)=M(i)*sind(beeta_5(i));
    Pr_1(i)=1+((2*gamma)/(gamma+1))*((M_n_1(i)*M_n_1(i))-1);
    Pr_2(i)=1+((2*gamma)/(gamma+1))*((M_n_2(i)*M_n_2(i))-1);
    Pr_3(i)=1+((2*gamma)/(gamma+1))*((M_n_3(i)*M_n_3(i))-1);
    Pr_4(i)=1+((2*gamma)/(gamma+1))*((M_n_4(i)*M_n_4(i))-1);
    Pr_5(i)=1+((2*gamma)/(gamma+1))*((M_n_5(i)*M_n_5(i))-1);
    C_p_1(i)=(2/(gamma*M_n_1(i)*M_n_1(i)))*(Pr_1(i)-1);
    C_p_2(i)=(2/(gamma*M_n_2(i)*M_n_2(i)))*(Pr_2(i)-1);
    C_p_3(i)=(2/(gamma*M_n_3(i)*M_n_3(i)))*(Pr_3(i)-1);
    C_p_4(i)=(2/(gamma*M_n_4(i)*M_n_4(i)))*(Pr_4(i)-1);
    C_p_5(i)=(2/(gamma*M_n_5(i)*M_n_5(i)))*(Pr_5(i)-1);
    x_1(i)=M(i)*theeta_1;
    y_1(i)=C_p_1(i)/(theeta_1*theeta_1);
    x_2(i)=M(i)*theeta_2;
    y_2(i)=C_p_2(i)/(theeta_2*theeta_2);
    x_3(i)=M(i)*theeta_3;
    y_3(i)=C_p_3(i)/(theeta_3*theeta_3);
    x_4(i)=M(i)*theeta_1;
    y_4(i)=C_p_4(i)/(theeta_4*theeta_4);
    x_5(i)=M(i)*theeta_1;
    y_5(i)=C_p_5(i)/(theeta_5*theeta_5)
end
figure(1)
plot(M,C_p_1,'r')
hold on;
plot(M,C_p_2,'g')
hold on;
plot(M,C_p_3,'b')
hold on;
plot(M,C_p_4,'k')
hold on;
plot(M,C_p_5,'y')
figure(2)
plot(x_1,y_1)
hold on;
plot(x_2,y_2)
plot(x_3,y_3)
plot(x_4,y_4)
plot(x_5,y_5)