clear; close all; clc;
load('given.mat')
A = Ebar\Abar
B = Ebar\Bbar
C =  [0 0 0 1;
      0 1 0 0]
wn_sp = 16;
dr_sp = .8;
n = -wn_sp*(dr_sp);
w = wn_sp*sqrt(1-(dr_sp)^2);
lam1 = n + 1i*w;
lam2 = conj(lam1);
lambda = [lam1, lam2];
[ol_vec, ol_lam] = eig(A);
% new_vec1 = [1+i*-1; -1+i*1; 0; 0];
% new_vec3 = [0; 0; 1+i*-1; -1+i*1];
% new_vec2 = conj(new_vec1);
% new_vec4 = conj(new_vec2);
% new_vec = [new_vec1, new_vec2, new_vec3, new_vec4];

% new_vec = [1-i, 1+i; 
%            -1+i, -1-i];
% new_vec = ol_vec([1,2],1:2)
% O = zeros(2,1); 
% D = [1 0 0 0;
%      0 1 0 0;];

% new_vec = ol_vec(:,1:2);
% O = zeros(4,1);
% D = eye(4);
% 
% I4 = eye(4);
% M1 = [(lambda(1)*I4-A), B; D, O];
% M1inv = pinv(M1);
% M2 = [(lambda(2)*I4-A), B; D, O];
% M2inv = pinv(M2);
% 
new_vec = [1-i, 1+i];
O = zeros(1,1); 
D = [1 0 0 0];

I4 = eye(4);
M1 = [(lambda(1)*I4-A), B; D, O];
M1inv = inv(M1);
M2 = [(lambda(2)*I4-A), B; D, O];
M2inv = inv(M2);

top_zeros = zeros(4,2);
des = [top_zeros; new_vec];
final1 = M1inv*des(:,1);
v1 = final1(1:4);
u1 = final1(5);
final2 = M2inv*des(:,2);
v2 = final2(1:4);
u2 = final2(5);
U = [u1, u2];
V = [v1, v2];
Kstruc = U*inv(C*V);
% Kstruc = -Kstruc;

% Kstruc = EigStructAssign(A, B, C, lambda(1:2), new_vec);
ssc = ss(A-B*Kstruc*C,B,C,0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
figure()
initial(ssc, [0, 0, 0, deg2rad(10)])
title("Eigenstructure Assignment gains Response to 10 Degree Initial Pitch")

vd1 = [eye(4), zeros(4,1)];
fprintf("Eigenstructure Assignment\n")
fprintf("Desired Short Period Poles:\n")
disp([lam1; lam2])
fprintf("Gains are K = [%f  %f] \n", Kstruc)
fprintf("Poles using Eigenstructure Assignment gains:\n")
[es_vec, es_lam] = eig(A-B*Kstruc*C);
disp(es_lam)
fprintf("Poles using Eigenstructure Assignment gains:\n")
disp(es_vec)