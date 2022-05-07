clear; close all; clc;

problem1()
% problem2()
problem3()

function problem1()
    A = [0,  1;
         0, -1];
    B = [0, 1]';
    Q = eye(2);
    R = 1;
    fprintf("Problem 1:\n")
    lqrprint(A, B, Q, R)
end

function problem2()
    A = [-1, 1;
          0, 2];
    B = [1, 0]';
    Q = eye(2);
    R = 1;
    fprintf("Problem 2:\n")
    lqrprint(A, B, Q, R)
end

function problem3() 
    A = [0, 1;
         0, 0];
    B = [0, 1]';
    Q = [1, 0;
         0, 1];
    R = 1;
    fprintf("Problem 3:\n")
    lqrprint(A, B, Q, R)
end

function lqrprint(A, B, Q, R)
    [K, S, CLP] = lqr(A, B, Q, R);
    fprintf("Gains:\n")
    disp(K)
    fprintf("P:\n")
    disp(S)
    fprintf("Closed Loop Poles:\n")
    disp(CLP)
end