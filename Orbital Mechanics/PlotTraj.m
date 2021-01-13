function PlotTraj(a,e,w,i,omega,theta,vtheta,fig)

if(nargin<8)
    fig=figure;
else
    figure(fig);
end

Re=6378.14; %earth radius

r=abs(a*(1-e^2)./(1+e*cos(vtheta))); %perifocal frame of the position vector
x=r.*cos(vtheta);
y=r.*sin(vtheta);

Q=perifocaltogeo(w,omega,i); %coordinate switch from perifocal to geocentric
rgeo=Q*[x;y;zeros(size(x))];
rx=rgeo(1,:);
ry=rgeo(2,:);
rz=rgeo(3,:);

%plotting
if(nargin<8)
plot3(rx,ry,rz,'-b')
grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
hold on, [x,y,z]=sphere(50);
surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')
else
plot3(rx,ry,rz,'-r')
grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
hold on, [x,y,z]=sphere(50);
surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')
end