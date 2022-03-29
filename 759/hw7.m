clear; close all; clc;

I = [.45, .75, .85, .65, .55];
epsilon = .2;
x = I;
fprintf("MaxNet\n")
for ii = 1:5
    [xnew, u] = maxnet(x, epsilon);
    fprintf("Time Step: %f\n", ii)
    fprintf("\tU:\n")
    disp(u)
    fprintf("\tY:\n")
    disp(xnew)
    x = xnew;
end

I = [.25, .85, .45, .55, .45];
R1 = 1;
w1 = .5;
R2 = 2;
w2 = -.3;
x = I;
fprintf("Mexican Hat\n")
for ii = 1:5
    [xnew, a] = mexican_hat(x, R1, w1, R2, w2);
    fprintf("Time Step: %f\n", ii)
    fprintf("\ta:\n")
    disp(a)
    fprintf("\tY:\n")
    disp(xnew)
    x = xnew;
end

w = [.3, .6, .1, .4, .8;
     .7, .9, .5, .3, .2];
alpha = .2;
x = [.5; .2];
fprintf("Kohonen\n")
for ii = 1:1
    [wnew, D] = kohonen(w, x, alpha);
    fprintf("Time Step: %f\n", ii)
    fprintf("\tD:\n")
    disp(D)
    fprintf("\tW:\n")
    disp(wnew)
    w = wnew;
end
