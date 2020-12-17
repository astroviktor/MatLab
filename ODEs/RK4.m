function y=RK4(t,x0,f)
%This function solves an ODE using the Runge-Kutta Fourth-Order Method
n=length(t);
dt=t(2)-t(1);
[r,~]=size(x0);
y=[x0,zeros(r,n-1)];
for i=1:n-1
    % Partial updates
    dy1=f(t(i),y(:,i));
    dy2=f(t(i)+dt/2,y(:,i)+dy1*dt/2);
    dy3=f(t(i)+dt/2,y(:,i)+dy2*dt/2);
    dy4=f(t(i)+dt,y(:,i)+dy3*dt);
    % Function update
    y(:,i+1)=y(:,i)+dt*(dy1+2*dy2+2*dy3+dy4)/6;
end
end