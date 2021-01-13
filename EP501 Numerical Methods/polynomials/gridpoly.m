function indexgrid=gridpoly(x,xprime)
%function that works over a grid. given a point and a grid,the function
%evaluates in which index of the grid that point falls
n=length(x);
for i=1:n-1
    if x(i)<=xprime && x(i+1)>=xprime
        index(i)=i;
    end
end
end
