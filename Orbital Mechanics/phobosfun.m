function val=phobosfun(x,y)
global a b lambda
val=((x-(1-lambda)).^2./a^2)+(y.^2./b^2)-1;
end