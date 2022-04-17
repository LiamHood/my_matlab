clear; close all; clc;

syms T alpha epsilon D L m g gamma lambda1 lambda2 V

f1 = T*cos(alpha + epsilon) - D - m*g*sin(gamma);
f2 = T*sin(alpha + epsilon) + L - m*g*cos(gamma);
L = V*sin(gamma);

H = L + lambda1*f1 + lambda2*f2;

