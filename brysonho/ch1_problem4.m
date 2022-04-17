clear; close all; clc;

syms lambda1
n = 2;
m = 2;
x = sym('x', [n 1]);
u = sym('u', [m 1]);
Q = sym('Q', [n n]);
R = sym('R', [m m]);
G = sym('G', [n m]);
c = sym('c', [n 1]);
L = (1/2)*x'*Q*x + (1/2)*u'*R*u;
f1 = x + G*u + c;
H = L + lambda1*f1;

dHdx = diff(H,x);
dHdu = diff(H,u);

sol = solve([dHdx; dHdu; f1]==[0; 0; 0], ...
    [x, u, lambda1]);
fprintf("solved for value\n")
disp("x = ")
disp(eval(sol.x))
disp("u = ")
disp(eval(sol.u))
disp("lambda1 = ")
disp(eval(sol.lambda1))

n = length(sol.x);
for ii = 1:n
    if sol.x(ii) > 0 && sol.y(ii) > 0 && sol.z(ii) > 0
        possible_maxval(ii,1) = eval(collect(subs(L,[x, y, z], ...
            [sol.x(ii), sol.y(ii), sol.z(ii)])));
    else
        possible_maxval(ii,1) = 0;
    end
    check(ii,1) = eval(subs(f1,[x, y, z], ...
        [sol.x(ii), sol.y(ii), sol.z(ii)]));
end
fprintf("Confirm restriction is 0\n")
disp(check)
[M,I] = max(possible_maxval);
fprintf("Max Volume: %f\n",M)
fprintf("Sides: x=%f  y=%f  z=%f\n", [sol.x(I), sol.y(I), sol.z(I)])