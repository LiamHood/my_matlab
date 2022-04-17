clear; close all; clc;

syms x u lambda
% a = 2;
% b = 3;
% c = 4;
% m = 5;

syms a b c m

L = (1/2)*(x^2/a^2 + u^2/b^2);
f = x + m*u - c ;
H = L + lambda*f;

dHdx = diff(H,x);
dHdu = diff(H,u);

sol = solve([dHdx; dHdu; f]==[0; 0; 0],[x,u,lambda]);
fprintf("solved for value\n")
disp(sol)

Lmin = collect(subs(L,[x, u, lambda], [sol.x, sol.u, sol.lambda]));
fprintf("Minimum value of L and cost\n")
disp(Lmin)