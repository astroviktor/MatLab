function finterp=interp2D(fvec,xvec,yvec,xprime,yprime)
    %matrix
    M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
    b=fvec(:);
    %gauss_elimination procedure
    [Mmod,ord]=Gauss_elim(M,b);
    aa=backsub(Mmod(ord,:));
    %interpolated value
   
    finterp=aa(1)+aa(2)*xprime+aa(3)*yprime+aa(4)*xprime*yprime;
        
end