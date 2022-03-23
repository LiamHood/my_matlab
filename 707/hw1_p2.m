fprintf("******Problem 2******\n")
% xdot = A*x + B*u
A = [8, -2, 3;
     0, 1, 4;
     7, 9, 10];
B = [2; 1; 5];
% y = C*x+D*u
C = [2, 1, 3];
D = 0;
disp("start")
% Transfer function G(s)
fprintf("Transfer Function\n")
syms s
G = collect(C*inv((s*eye(3)-A))*B);
fprintf("From formula\n")
fprintf("G(s) = %s\n", G)
fprintf("Built in Function\n")
[b,a] = ss2tf(A,B,C,D);
G_builtin = collect((b(1)*s^3 + b(2)*s^2 + b(3)*s + b(4))/...
    (a(1)*s^3 + a(2)*s^2 + a(3)*s + a(4)));
fprintf("G(s) = %s\n", G_builtin)
% Poles
fprintf("\nPoles\n")
poles_builtin = eig(A);
poles = vpasolve(det(s*eye(3)-A)==0,s);
fprintf("From formula\n")
fprintf("Poles %f, %f, %f \n", poles)
fprintf("Built in Function\n")
fprintf("Poles %f, %f, %f \n", poles_builtin)