clear; close all; clc;

syms x y lambda
% syms a b
a = 2;
b = 3;

L = 4*(x + y);
f = x^2/a^2 + y^2/b^2 - 1;
H = L + lambda*f;

dHdx = diff(H,x);
dHdy = diff(H,y);

sol = solve([dHdx; dHdy; f]==[0; 0; 0], ...
    [x, y, lambda]);
fprintf("solved for value\n")
disp("x = ")
disp(eval(sol.x))
disp("y = ")
disp(eval(sol.y))
disp("lambda = ")
disp(eval(sol.lambda))

n = length(sol.x);
for ii = 1:n
    if sol.x(ii) > 0 && sol.y(ii) > 0
        possible_maxval(ii,1) = eval(collect(subs(L,[x, y], ...
            [sol.x(ii), sol.y(ii)])));
    else
        possible_maxval(ii,1) = 0;
    end
    check(ii,1) = eval(subs(f,[x, y], ...
        [sol.x(ii), sol.y(ii)]));
end
fprintf("Confirm restriction is 0\n")
disp(check)
[M,I] = max(possible_maxval);
fprintf("Max Perimeter: %f\n",M)
fprintf("Sides: x=%f  y=%f\n", [sol.x(I), sol.y(I)])