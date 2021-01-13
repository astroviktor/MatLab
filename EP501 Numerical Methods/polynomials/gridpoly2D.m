function [indexgridx,indexgridy]=gridpoly2D(x,xprime,y,yprime)
%function that works over a grid. given a point and a grid,the function
%evaluates in which index of the grid that point falls
nx=length(x);
for i=2:nx
    if x(i-1)<=xprime && x(i)>=xprime
        indexgridx=i-1;
        break
    end
end
ny=length(y);
for i=2:ny
    if y(i-1)<=yprime && y(i)>=yprime
        indexgridy=i-1;
        break
    end
end
end
