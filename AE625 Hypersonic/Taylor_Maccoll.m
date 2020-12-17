function y=Taylor_Maccoll(f,x0,y0,z0,maxit,dx)

%{ 
This function solves taylor_maccoll equations using runge-kutta numerical
method

-f           function to handle
-x0,y0,z0    initial conditions
-dx          infinitesimal increment
-maxit       maximum number of iterations (to avoid loops)
%}

i=1;

while i<=maxit
    
    m1=dx*f(x0,y0,z0);
    k1=dx*z0;
    m2=dx*f(x0+dx/2,y0+k1/2,z0+m1/2);
    k2=dx*z0;
    m3=dx*f(x0+dx/2,y0+k2/2,z0+m2/2);
    k3=dx*z0;
    m4=dx*f(x0+dx,y0+k3,z0+m3);
    
    
    
    
    

