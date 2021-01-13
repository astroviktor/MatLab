function [x,it]=secanti(f,x0,x1,toll,maxit)

arresto=1;
it=0;

while arresto>toll && abs(f(x1))>=toll && it<=maxit
    c=x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    arresto=abs(c-x1)/abs(c);
    
    x0=x1;
    x1=c;
    
   it=it+1;
end
x=c;
it=it-1;
