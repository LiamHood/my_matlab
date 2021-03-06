clear; close all; clc;
load('given.mat')
A = [-1.341, .9933, 0 , -.1689, -.2518;
     43.223, -.8693, 0, -17.251, -1.5766;
     1.34, .0067, 0, .1689, .2518;
     0, 0, 0, -20, 0;
     0, 0, 0, 0, -20];
B = [0, 0, 0, 20, 0;
     0, 0, 0, 0, 20]';
C = [ 0 1 0 0 0;
     47.76, -.268, 0, -4.56, 4.45;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
lam1 = -5.6 + 1i*4.2;
lam2 = conj(lam1);
lambda = [lam1, lam2, -1, -19, -19.5];
[ol_vec, ol_lam] = eig(A);
% new_vec1 = [1+i*-1; -1+i*1; 0; 0];
% new_vec3 = [0; 0; 1+i*-1; -1+i*1];
% new_vec2 = conj(new_vec1);
% new_vec4 = conj(new_vec2);
% new_vec = [new_vec1, new_vec2, new_vec3, new_vec4];

new_vec = [1-i, 1+i; 
           -1+i, -1-i;
           0, 0];
O = zeros(3,2); 
D = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0];
I = eye(5);
M1 = [(lambda(1)*I-A), B; D, O];
M1inv = pinv(M1);
M2 = [(lambda(2)*I-A), B; D, O];
M2inv = pinv(M2);

new_vec3 = [0, 1]';
O = zeros(2,2); 
D = [0 1 0 0 0;
     0 0 1 0 0];
M3 = [(lambda(3)*I-A), B; D, O];
M3inv = pinv(M3);

new_vec4 = [1];
O = zeros(1,2); 
D = [0 0 0 1 0];
M4 = [(lambda(4)*I-A), B; D, O];
M4inv = pinv(M4);

new_vec5 = [1];
O = zeros(1,2); 
D = [0 0 0 0 1];
M5 = [(lambda(5)*I-A), B; D, O];
M5inv = pinv(M5);

top_zeros = zeros(5,1);
des1 = [top_zeros; new_vec(:,1)];
des2 = [top_zeros; new_vec(:,2)];
des3 = [top_zeros; new_vec3];
des4 = [top_zeros; new_vec4];
des5 = [top_zeros; new_vec5];
final1 = M1inv*des1;
v1 = final1(1:5);
u1 = final1(6:7);
final2 = M2inv*des2;
v2 = final2(1:5);
u2 = final2(6:7);
final3 = M3inv*des3;
v3 = final3(1:5);
u3 = final3(6:7);
final4 = M4inv*des4;
v4 = final4(1:5);
u4 = final4(6:7);
final5 = M5inv*des5;
v5 = final5(1:5);
u5 = final5(6:7);
U = [u1, u2];
V = [v1, v2];
Kstruc = U*pinv(C*V);
% Kstruc = -Kstruc;

% Kstruc = EigStructAssign(A, B, C, lambda(1:2), new_vec);
ssc = ss(A-B*Kstruc*C,zeros(5,1),C,0);
figure()
initial(ssc, [0, 0, 0, 0, deg2rad(10)])
title("Eigenstructure Assignment gains Response to 10 Degree Initial Pitch")

vd1 = [eye(4), zeros(4,1)];
fprintf("Eigenstructure Assignment\n")
fprintf("Desired Short Period Poles:\n")
disp([lam1; lam2])
% fprintf("Gains are K = [%f  %f] \n", Kstruc)
disp(Kstruc)
fprintf("Poles using Eigenstructure Assignment gains:\n")
[es_vec, es_lam] = eig(A-B*Kstruc*C);
disp(es_lam)
fprintf("vectors using Eigenstructure Assignment gains:\n")
disp(es_vec)