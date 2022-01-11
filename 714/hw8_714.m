clear; close all; clc;
A = [0, 1, 0, 2;
    5, 5, 5, 5;
    20, 10, 16, 8;
    51, 17, 39, 13];
b = [3; 17; 42; 93];

answer = inv(A)*b