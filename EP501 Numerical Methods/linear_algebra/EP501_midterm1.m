%% EP 501 - Numerical Methods
%% Midterm
% Vittorio Baraldi

clear all; close all; clc

%% Exercise 1
% Part a)
load C:\Users\Vittorio\EP501_assignments\assignments\HW2\iterative_testproblem.mat
%solving tridiagonal system using thomas algorithm
x=tridiag(Ait,bit);
%verifying solution
xcheck=Ait\bit;
%outputs
disp('Thomas algorithm solution:')
disp(x)
disp('MatLab built-in function solution:')
disp(xcheck)