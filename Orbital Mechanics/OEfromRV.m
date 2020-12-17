function [a,e,i,omega,theta,w]=OEfromRV(r,v,mu)

toll=2*eps;

rm=sqrt(r'*r);
vm=sqrt(v'*v);
h=cross(r,v); %angular momentum
hm=sqrt(h'*h);
i=acos(h(3)/hm); %inclination
N=cross([0;0;1],h);
Nm=sqrt(N'*N);
%checking particular cases for i=0 and i=pi
if Nm<toll
    N=[1;0;0];
    Nm=1.0;
end

%right ascension
if (abs(i)<toll || abs(i-pi)<toll)
    omega=0;
elseif N(2)>=2
    omega=acos(N(1)/Nm);
else
    omega=2*pi-acos(N(1)/Nm);
end

%eccentricity
vr=v'*r/rm;
ev=1/mu*((vm^2-mu/rm)*r-vr*rm*v);
e=sqrt(ev'*ev);

%semimayor axis
a=hm^2/mu/(1-e^2);

%setting parameters for special cases
equatorial=(abs(i)<toll || abs(i-pi)<toll);
circular=e<toll;

if (circular) %argument of perigee
    w=0;
else
    if (equatorial)
    arg=(N'*ev)/(e*Nm);
    if(arg>=1)
        w=0;
    elseif(arg<=-1);
        w=pi;
    else
        w=acos(arg);
    end
    %retrograde case:
    if(abs(i)<toll && ev(2)<0)
      w = 2*pi-w;
    elseif(abs(i-pi)<toll && ev(2)>0)
      w = 2*pi-w;
    end
else
    arg=(N'*ev)/(e*Nm);
    if(arg>=1)
        w=0;
    elseif(arg<=-1);
        w=pi;
    else
        w=acos(arg);
    end
    %retrograde case:
    if(ev(3)<0)
      w = 2*pi-w;
    elseif(abs(i-pi)<toll && ev(2)>0)
      w = 2*pi-w;
    end
  end
end

%True anomaly
if (circular && equatorial)
    theta=acos(r(1)/rm);
    if r(2)<0
        theta=2*pi-theta;
    end
else
    if (circular)
        eccdirection=N/Nm;
        arg=(eccdirection'*r)/rm;
    else
        eccdirection=ev/e;
        arg=(eccdirection'*r)/rm;
    end
if(arg>=1.0)
    arg=1.0;
  elseif(arg<=-1.0)
    arg=-1.0;
end
if(vr>=0)
    theta=acos(arg); 
else
    theta=2*pi-acos(arg);
end
end