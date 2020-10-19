%% EP 501
%% Vittorio Baraldi
% Homework 1
clear all
close all
%% Ex. 3

load C:\Users\Vittorio\EP501_matlab\linear_algebra\testproblem.mat

[matrix,ord,d]=Gausselim_det(A,b);

dcheck=det(A);

disp('Determinant obtained from modified Gauss elimination function: ')
disp(d)
disp('Determinant obtained from built-in MatLab function: ')
disp(dcheck)