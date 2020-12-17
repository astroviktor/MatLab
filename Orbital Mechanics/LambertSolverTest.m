function LambertSolverTest(r1,v1,r2,v2,mu)

[a,e,i,omega,theta,w]=OEfromRV(r1,v1,mu);
[a2,e2,i2,omega2,theta2,w2]=OEfromRV(r2,v2,mu);

if(norm([1,cos([i,omega,w]),e]-[a2/a,cos([i2,omega2,w2]),e2])<10e-5)
    fprintf("Your lambert solver works well!\n\n");
else
    frpintf("Orbital elements solutions DON'T match! Review your calculations\n\n");
end