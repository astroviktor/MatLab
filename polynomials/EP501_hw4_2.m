%% EP 501 - Project 4
% Vittorio Baraldi
clear all; close all; clc

%% Exercise 2

%% Part a) thru c)

load test_interp.mat
addpath ../linear_algebra

n=length(xgi);    sol=[]; xf=[]; yf=[]; fvec=[];
for i=1:n
    %evaluating indexes such that the point to which data are to be
    %interpolated 
    [indexgridx(i),indexgridy(i)]=gridpoly2D(xg,xgi(i),yg,ygi(i));
    %4 points in which data are to be interpolated
    xf=[xf;xgi(indexgridx(i)) xgi(indexgridx(i)+1)];
    yf=[yf;ygi(indexgridy(i)) ygi(indexgridy(i)+1)];   
    %grid mesh
    [x,y]=meshgrid(xf(i,:),yf(i,:));
    %vectorizing elements
    xvec=x(:);
    yvec=y(:);
    %picking the function points 
    fvec=[f2D(indexgridx(i)),f2D(indexgridx(i)+1);f2D(indexgridy(i)),f2D(indexgridy(i)+1)];
    ff=fvec(:);
    %interpolating
    fun=interp2D(ff,xvec,yvec,xgi(i),ygi(i));
    sol=[sol;fun];
end
%matlab built-in function
finterp=interp2(xg,yg,f2D,xgi,ygi);
%outputs and plots
figure(1)
plot(xgi,finterp)
title('Built-in MatLab interpolation')
xlabel('x')
ylabel('f(x)')
figure(2)
plot(xgi,sol)
title('User interpolation')
xlabel('x')
ylabel('f(x)')