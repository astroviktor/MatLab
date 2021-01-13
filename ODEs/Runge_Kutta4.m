function y=Runge_Kutta4(y0,alpha,t)
lt=numel(t);
dt=t(2)-t(1);
y=zeros(1,lt);
y(1)=y0;
if length(alpha)>1
        for n=1:lt-1
  
            dy1=-alpha(n)*dt*y(n);
            dy2=-(alpha(n)*dt)*y(n)*(1-(alpha(n)*dt)/2);
            dy3=-(alpha(n)*dt)*y(n)*(1-(alpha(n)*dt)/2+(alpha(n)*dt)^2/4);
            dy4=-(alpha(n)*dt)*y(n)*(1-(alpha(n)*dt)+(alpha(n)*dt)^2/2-(alpha(n)*dt)^3/4);
            y(n+1)=y(n)+1/6*(dy1+2*dy2+2*dy3+dy4);
        end
        
else
        for n=1:lt-1
            dy1=-alpha*dt*y(n);
            dy2=-(alpha*dt)*y(n)*(1-(alpha*dt)/2);
            dy3=-(alpha*dt)*y(n)*(1-(alpha*dt)/2+(alpha*dt)^2/4);
            dy4=-(alpha*dt)*y(n)*(1-(alpha*dt)+(alpha*dt)^2/2-(alpha*dt)^3/4);
            y(n+1)=y(n)+1/6*(dy1+2*dy2+2*dy3+dy4);
        end
end

end