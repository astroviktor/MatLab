S=47;
T=40298*2;
c0=0.032;
W=103047;
e=0.87;
AR=6.5;
p=1.225;
k=1/(pi*e*AR);
for v=0:10:500
    
    cl=(2*W)/(p*v^2*S);
    cd=c0+k*cl^2;
    T=W/(cl/cd);
    P=T*v;
    plot(v,P,'-r');
end

