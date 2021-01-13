%% EP 501
%% Vittorio Baraldi
% Homework 1
clear all
close all

%% Ex. 2
load C:\Users\Vittorio\EP501_matlab\linear_algebra\testproblem.mat
%finding the inverse matrix of testproblem.mat using the implemented gauss
%jordan elimination function
c=eye(size(A,1));
x=gauss_jordan_elim(A,c);
%checking results with the built-in MatLab function
xcheck=inv(A);

disp('Solution for Gauss-Jordan elimination: ')
disp(x)
disp('Solution for built-in MatLab function: ')
disp(xcheck)