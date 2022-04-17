clear; close all; clc;

syms x y z lambda1 lambda2
% syms a b c
a = 2;
b = 3;
c = 5;

L = 8*x*y*z;
f1 = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1;
H = L + lambda1*f1;

dHdx = diff(H,x);
dHdy = diff(H,y);
dHdz = diff(H,z);

sol = solve([dHdx; dHdy; dHdz; f1]==[0; 0; 0; 0], ...
    [x, y, z, lambda1]);
fprintf("solved for value\n")
disp("x = ")
disp(eval(sol.x))
disp("y = ")
disp(eval(sol.y))
disp("z = ")
disp(eval(sol.z))
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