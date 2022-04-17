clear; close all; clc;

syms x y z lambda1 lambda2

L = x^2 + y^2 + z^2;
f1 = x + 2*y + 3*z - 10;
f2 = x - y + 2*z - 1;
H = L + lambda1*f1 + lambda2*f2;

dHdx = diff(H,x);
dHdy = diff(H,y);
dHdz = diff(H,z);

sol = solve([dHdx; dHdy; dHdz; f1; f2]==[0; 0; 0; 0; 0], ...
    [x, y, z, lambda1, lambda2]);
fprintf("solved for value\n")
disp(sol)

Lmin = collect(subs(L,[x, y, z, lambda1, lambda2], ...
    [sol.x, sol.y, sol.z, sol.lambda1, sol.lambda2]));
fprintf("Minimum value of L and cost\n")
disp(Lmin)