rc = [  5.79e7; %mercury
        1.08e8; %venus
        149.6e6; %earth 
        2.28e8; %mars
        7.78e8; %jupiter
        1.43e9; %saturn
        2.87e9; %uranus
        4.5e9];%neptune 
re=149.6e6;
mus=1.327e11;
n=OrbRate(rc,mus);
ne=n(3);

    Tsyn=(2*pi)./abs(n'-ne)/86400;

for i=1:length(rc)
    k(i)=HohmannDV(re,rc(i),mus)
end
k

Tsyn
