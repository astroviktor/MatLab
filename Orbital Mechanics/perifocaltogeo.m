function Q=perifocaltogeo(w,o,i)

q1=[cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];
q2=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
q3=[cos(o) sin(o) 0; -sin(o) cos(o) 0; 0 0 1];

Q=q1*q2*q3;
Q=Q';







