clc;
clear all;
close all;
gamma=1.4;
M=1.5:0.1:7 
for i=1:length(M)
    theta_1=5;
    theta_2=10;
    theta_3=15;
    theta_4=20;
    theta_5=30;
    C_p_1(i)= 2*theta_1*theta_1*[((gamma+1)/4) + sqrt((((gamma+1)/4)^2)+(1/((theta_1*M(i))^2)))];
    C_p_2(i)= 2*theta_2*theta_2*[((gamma+1)/4) + sqrt((((gamma+1)/4)^2)+(1/((theta_2*M(i))^2)))];
    C_p_3(i)= 2*theta_3*theta_3*[((gamma+1)/4) + sqrt((((gamma+1)/4)^2)+(1/((theta_3*M(i))^2)))];
    C_p_4(i)= 2*theta_4*theta_4*[((gamma+1)/4) + sqrt((((gamma+1)/4)^2)+(1/((theta_4*M(i))^2)))];
    C_p_5(i)= 2*theta_5*theta_5*[((gamma+1)/4) + sqrt((((gamma+1)/4)^2)+(1/((theta_5*M(i))^2)))];
    y_1(i)=C_p_1(i)/(theta_1*theta_1);
    y_2(i)=C_p_2(i)/(theta_2*theta_2);
    y_3(i)=C_p_3(i)/(theta_3*theta_3);
    y_4(i)=C_p_4(i)/(theta_4*theta_4);
    y_5(i)=C_p_5(i)/(theta_5*theta_5);
    x_1(i)=M(i)*theta_1;
    x_2(i)=M(i)*theta_2;
    x_3(i)=M(i)*theta_3;
    x_4(i)=M(i)*theta_4;
    x_5(i)=M(i)*theta_5;
end
plot(x_1,y_1)
hold on;
plot(x_2,y_2)
plot(x_3,y_3)
plot(x_4,y_4)
plot(x_5,y_5)