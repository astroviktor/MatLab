function F=state(t,U)

q=1.6e-19;
m=1.67e-27;
B=50000e-9;
F=[U(3);U(4);q*(U(4)*B*(1+U(2)/2))/m;-q*(U(3)*B*(1+U(2)/2))/m];
end