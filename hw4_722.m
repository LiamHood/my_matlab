%% Problem 1
clear; close all; clc;
syms x1 x h del1 del2 del3 del4 del5 del6 phi;
xmatrix = [ 1, x1, x1^2, x1^3, x1^4, x1^5;
            0, 1, 2*x1, 3*x1^2, 4*x1^3, 5*x1^4;
            0, 0, 2, 6*x1, 12*x1^2, 20*x1^3;
            1, x1+h, (x1+h)^2, (x1+h)^3, (x1+h)^4, (x1+h)^5;
            0, 1, 2*(x1+h), 3*(x1+h)^2, 4*(x1+h)^3, 5*(x1+h)^4;
            0, 0, 2, 6*(x1+h), 12*(x1+h)^2, 20*(x1+h)^3];
delvec = [del1;del2;del3;del4;del5;del6];
cvec = inv(xmatrix)*delvec;
wcluttered = cvec(1) + cvec(2)*x + cvec(3)*x^2 + cvec(4)*x^3 + cvec(5)*x^4 + cvec(6)*x^5;
w = collect(wcluttered, [del1, del2, del3, del4, del5, del6]);
phi1 = simplify(solve( del1*phi == subs(w,[del2, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi2 = simplify(solve( del2*phi == subs(w,[del1, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi3 = simplify(solve( del3*phi == subs(w,[del1, del2, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi4 = simplify(solve( del4*phi == subs(w,[del1, del2, del3, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi5 = simplify(solve( del5*phi == subs(w,[del1, del2, del3, del4, del6, x1], [0 0 0 0 0 0]),phi));
phi6 = simplify(solve( del6*phi == subs(w,[del1, del2, del3, del4, del5, x1], [0 0 0 0 0 0]),phi));
fprintf("Problem 1 \n")
fprintf("phi1 = %s \n",phi1)
fprintf("phi2 = %s \n",phi2)
fprintf("phi3 = %s \n",phi3)
fprintf("phi4 = %s \n",phi4)
fprintf("phi5 = %s \n",phi5)
fprintf("phi6 = %s \n",phi6)
fprintf("phi1, phi2, phi3 are for node 1\n")

%% Problem 2
clear;
syms x1 x h del1 del2 del3 del4 del5 del6 phi;
xmatrix = [ 1, x1, x1^2, x1^3, x1^4, x1^5;
            0, -1, -2*x1, -3*x1^2, -4*x1^3, -5*x1^4;
            1, (x1+h/2), (x1+h/2)^2, (x1+h/2)^3, (x1+h/2)^4, (x1+h/2)^5;
            0, -1, -2*(x1+h/2), -3*(x1+h/2)^2, -4*(x1+h/2)^3, -5*(x1+h/2)^4;
            1, (x1+h), (x1+h)^2, (x1+h)^3, (x1+h)^4, (x1+h)^5;
            0, -1, -2*(x1+h), -3*(x1+h)^2, -4*(x1+h)^3, -5*(x1+h)^4];
delvec = [del1;del2;del3;del4;del5;del6];
cvec = inv(xmatrix)*delvec;
wcluttered = cvec(1) + cvec(2)*x + cvec(3)*x^2 + cvec(4)*x^3 + cvec(5)*x^4 + cvec(6)*x^5;
w = collect(wcluttered, [del1, del2, del3, del4, del5, del6]);
phi1 = simplify(solve( del1*phi == subs(w,[del2, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi2 = simplify(solve( del2*phi == subs(w,[del1, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi3 = simplify(solve( del3*phi == subs(w,[del1, del2, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi4 = simplify(solve( del4*phi == subs(w,[del1, del2, del3, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi5 = simplify(solve( del5*phi == subs(w,[del1, del2, del3, del4, del6, x1], [0 0 0 0 0 0]),phi));
phi6 = simplify(solve( del6*phi == subs(w,[del1, del2, del3, del4, del5, x1], [0 0 0 0 0 0]),phi));
fprintf("\nProblem 2 \n")
fprintf("phi1 = %s \n",phi1)
fprintf("phi2 = %s \n",phi2)
fprintf("phi3 = %s \n",phi3)
fprintf("phi4 = %s \n",phi4)
fprintf("phi5 = %s \n",phi5)
fprintf("phi6 = %s \n",phi6)
fprintf("phi1, phi2 are for node 1\n")


%% Problem 3
clear;
syms x1 x h del1 del2 del3 del4 del5 del6 phi q1 q2;
xmatrix = [ 1, x1, x1^2, x1^3, x1^4, x1^5;
            0, 1, 2*x1, 3*x1^2, 4*x1^3, 5*x1^4;
            0, 0, 2, 6*x1, 12*x1^2, 20*x1^3;
            1, x1+h, (x1+h)^2, (x1+h)^3, (x1+h)^4, (x1+h)^5;
            0, 1, 2*(x1+h), 3*(x1+h)^2, 4*(x1+h)^3, 5*(x1+h)^4;
            0, 0, 2, 6*(x1+h), 12*(x1+h)^2, 20*(x1+h)^3];
delvec = [del1;del2;del3;del4;del5;del6];
cvec = inv(xmatrix)*delvec;
wcluttered = cvec(1) + cvec(2)*x + cvec(3)*x^2 + cvec(4)*x^3 + cvec(5)*x^4 + cvec(6)*x^5;
w = collect(wcluttered, [del1, del2, del3, del4, del5, del6]);
phi1 = simplify(solve( del1*phi == subs(w,[del2, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi2 = simplify(solve( del2*phi == subs(w,[del1, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi3 = simplify(solve( del3*phi == subs(w,[del1, del2, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi4 = simplify(solve( del4*phi == subs(w,[del1, del2, del3, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi5 = simplify(solve( del5*phi == subs(w,[del1, del2, del3, del4, del6, x1], [0 0 0 0 0 0]),phi));
phi6 = simplify(solve( del6*phi == subs(w,[del1, del2, del3, del4, del5, x1], [0 0 0 0 0 0]),phi));
fprintf("\nProblem 3 \n")
fprintf("phi1 = %s \n",phi1)
fprintf("phi2 = %s \n",phi2)
fprintf("phi3 = %s \n",phi3)
fprintf("phi4 = %s \n",phi4)
fprintf("phi5 = %s \n",phi5)
fprintf("phi6 = %s \n",phi6)
q = [   int(((q1+x*(q2-q1)/h)*phi1),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi2),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi3),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi4),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi5),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi6),x,0,h)];
for ii = 1:6
    fprintf("q_%i = %s \n",ii, q(ii))
end
fprintf("q_1, q_2, q_3 are the components of the source vector for node 1\n")

%% Problem 4
clear;
syms x1 x h del1 del2 del3 del4 del5 del6 phi q1 q2;
xmatrix = [ 1, x1, x1^2, x1^3, x1^4, x1^5;
            0, -1, -2*x1, -3*x1^2, -4*x1^3, -5*x1^4;
            1, (x1+h/2), (x1+h/2)^2, (x1+h/2)^3, (x1+h/2)^4, (x1+h/2)^5;
            0, -1, -2*(x1+h/2), -3*(x1+h/2)^2, -4*(x1+h/2)^3, -5*(x1+h/2)^4;
            1, (x1+h), (x1+h)^2, (x1+h)^3, (x1+h)^4, (x1+h)^5;
            0, -1, -2*(x1+h), -3*(x1+h)^2, -4*(x1+h)^3, -5*(x1+h)^4];
delvec = [del1;del2;del3;del4;del5;del6];
cvec = inv(xmatrix)*delvec;
wcluttered = cvec(1) + cvec(2)*x + cvec(3)*x^2 + cvec(4)*x^3 + cvec(5)*x^4 + cvec(6)*x^5;
w = collect(wcluttered, [del1, del2, del3, del4, del5, del6]);
phi1 = simplify(solve( del1*phi == subs(w,[del2, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi2 = simplify(solve( del2*phi == subs(w,[del1, del3, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi3 = simplify(solve( del3*phi == subs(w,[del1, del2, del4, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi4 = simplify(solve( del4*phi == subs(w,[del1, del2, del3, del5, del6, x1], [0 0 0 0 0 0]),phi));
phi5 = simplify(solve( del5*phi == subs(w,[del1, del2, del3, del4, del6, x1], [0 0 0 0 0 0]),phi));
phi6 = simplify(solve( del6*phi == subs(w,[del1, del2, del3, del4, del5, x1], [0 0 0 0 0 0]),phi));
fprintf("\nProblem 4 \n")
fprintf("phi1 = %s \n",phi1)
fprintf("phi2 = %s \n",phi2)
fprintf("phi3 = %s \n",phi3)
fprintf("phi4 = %s \n",phi4)
fprintf("phi5 = %s \n",phi5)
fprintf("phi6 = %s \n",phi6)
q = [   int(((q1+x*(q2-q1)/h)*phi1),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi2),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi3),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi4),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi5),x,0,h);
        int(((q1+x*(q2-q1)/h)*phi6),x,0,h)];
for ii = 1:6
    fprintf("q_%i = %s \n",ii, q(ii))
end
fprintf("The source vector components for node 1 are q_1 and q_2\n")
fprintf("The source vector components for node 2 are q_3 and q_4\n")
fprintf("The source vector components for node 1 are q_5 and q_6\n")

%% Problem 5
clear;
syms xbar x
L = 144;
h = L/2;
E = 30e6;
b = 3;
c = 2;
I = b*(2*c)^3/12;
EI = E*I;
q0 = 50;

phi(1) = 1 - 3.*(xbar/h).^2 + 2.*(xbar/h).^2;
phi(2) = -xbar.*(1-xbar/h).^2;
phi(3) = 3.*(xbar/h).^2 - 2.*(xbar/h).^3;
phi(4) = -xbar.*((xbar/h).^2 - (xbar/h));

q1 = zeros(4,1);
q2 = zeros(4,1);
for ii = 1:4
    q1(ii) = int(q0*phi(ii),xbar,0,h);
    q2(ii) = int(q0*phi(ii),xbar,0,h);
end
q = [q1(1); q1(2); q1(3) + q2(1); q1(4) + q2(2); q2(3); q2(4)];

k = @(ii,jj) int(EI*diff(diff(phi(ii)))*diff(diff(phi(jj))),xbar,0,h);
K1 = zeros(4);
K2 = zeros(4);
for ii = 1:4
    for jj = 1:4
        K1(ii,jj) = k(ii,jj);
        K2(ii,jj) = k(ii,jj);
    end
end
Kg = zeros(6);
for ii = 1:4
    for jj = 1:4
        Kg(ii,jj) = K1(ii,jj);
    end
end

for ii = 1:4
    for jj = 1:4
        Kg(ii+2,jj+2) = Kg(ii+2,jj+2) + K2(ii,jj);
    end
end

Qtop = Kg(1:3,1:3)*[0;0;0]-q(1:3);
Q = [Qtop; 0; 0; 0];
DELbot = inv(Kg(4:6,4:6))*(q(1:3)+Q(4:6));
DEL = [0;0;0;DELbot];

w1 = 0;
w2 = 0;
for ii = 1:4
    w1 = w1 + DEL(ii).*phi(ii);
    w2 = w2 + DEL(ii+2).*phi(ii);
end
w1 = simplify(w1);
w2 = simplify(w2);
M1 = EI*diff(diff(w1));
M2 = EI*diff(diff(w2));
sigma1 = M1*c/I;
sigma2 = M2*c/I;
fprintf("\nProblem5\n")
fprintf("Q_equil force at left: %f lb\n", Q(1))
fprintf("Q_equil moment at left: %f in*lb\n", Q(2))
fprintf("Q_equil force at right: %f lb\n", Q(5))
fprintf("Q_equil moment at right: %f in*lb\n", Q(6))
Qdef(1) = subs(diff(M1),xbar,0);
Qdef(2) = subs(M1,xbar,0);
Qdef(3) = -subs(diff(M2),xbar,L/2);
Qdef(4) = -subs(M2,xbar,L/2);
fprintf("Q_def force at left: %f lb\n", Qdef(1))
fprintf("Q_def moment at left: %f in*lb\n", Qdef(2))
fprintf("Q_def force at right: %f lb\n", Qdef(3))
fprintf("Q_def moment at right: %f in*lb\n", Qdef(4))

xvec1 = linspace(0,L/2);
xvec2 = linspace(L/2,L);
figure
hold on
title("Problem 5 w(x)")
plot(xvec1, subs(w1, xbar, xvec1))
plot(xvec2, subs(w2, xbar, xvec2-L/2))
xlabel("x (in)")
ylabel("Vertical displacement (in)")
hold off
figure
hold on
title("Problem 5 \sigma(x)")
plot(xvec1, subs(sigma1, xbar, xvec1))
plot(xvec2, subs(sigma2, xbar, xvec2-L/2))
xlabel("x (in)")
ylabel("Bending Stress (psi)")
hold off

fprintf("vertical dispalacement for 0<=x<=L/2 is w(x)=%s\n", subs(w1,xbar,x))
fprintf("vertical dispalacement for L/2<=x<=L is w(x)=%s\n", simplify(subs(w2,xbar,x-L/2)))
fprintf("bending stress for 0<=x<=L/2 is sigma(x)=%s\n", subs(sigma1,xbar,x))
fprintf("bending stress for L/2<=x<=L is sigma(x)=%s\n", simplify(subs(sigma2,xbar,x-L/2)))
%% Problem 6
clear;
syms xbar x
L = 144;
h = L/2;
E = 30e6;
b = 3;
c = 2;
I = b*(2*c)^3/12;
EI = E*I;

q0 = 50;
Ks = 2000;
M0 = 10000;
F0 = 2000;

phi(1) = 1 - 3.*(xbar/h).^2 + 2.*(xbar/h).^2;
phi(2) = -xbar.*(1-xbar/h).^2;
phi(3) = 3.*(xbar/h).^2 - 2.*(xbar/h).^3;
phi(4) = -xbar.*((xbar/h).^2 - (xbar/h));

q1 = zeros(4,1);
q2 = zeros(4,1);
for ii = 1:4
    q1(ii) = int(q0*phi(ii),xbar,0,h);
    q2(ii) = 0;
end
q = [q1(1); q1(2); q1(3) + q2(1); q1(4) + q2(2); q2(3); q2(4)];

k = @(ii,jj) int(EI*diff(diff(phi(ii)))*diff(diff(phi(jj))),xbar,0,h);
K1 = zeros(4);
K2 = zeros(4);
for ii = 1:4
    for jj = 1:4
        K1(ii,jj) = k(ii,jj);
        K2(ii,jj) = k(ii,jj);
    end
end
K2(3,3) = K2(3,3) + Ks;
Kg = zeros(6);
for ii = 1:4
    for jj = 1:4
        Kg(ii,jj) = K1(ii,jj);
    end
end

for ii = 1:4
    for jj = 1:4
        Kg(ii+2,jj+2) = Kg(ii+2,jj+2) + K2(ii,jj);
    end
end

Qtop = Kg(1:3,1:3)*[0;0;0]-q(1:3);
Q = [Qtop; -M0; F0; 0];
DELbot = inv(Kg(4:6,4:6))*(q(1:3)+Q(4:6));
DEL = [0;0;0;DELbot];
Q(5) = Q(5) - Ks*DEL(5);

w1 = 0;
w2 = 0;
for ii = 1:4
    w1 = w1 + DEL(ii).*phi(ii);
    w2 = w2 + DEL(ii+2).*phi(ii);
end
w1 = simplify(w1);
w2 = simplify(w2);
M1 = EI*diff(diff(w1));
M2 = EI*diff(diff(w2));
sigma1 = M1*c/I;
sigma2 = M2*c/I;
fprintf("\nProblem 6\n")
fprintf("Q_equil force at left: %f lb\n", Q(1))
fprintf("Q_equil moment at left: %f in*lb\n", Q(2))
fprintf("Q_equil force at right: %f lb\n", Q(5))
fprintf("Q_equil moment at right: %f in*lb\n", Q(6))
Qdef(1) = subs(diff(M1),xbar,0);
Qdef(2) = subs(M1,xbar,0);
Qdef(3) = -subs(diff(M2),xbar,L/2);
Qdef(4) = -subs(M2,xbar,L/2);
fprintf("Q_def force at left: %f lb\n", Qdef(1))
fprintf("Q_def moment at left: %f in*lb\n", Qdef(2))
fprintf("Q_def force at right: %f lb\n", Qdef(3))
fprintf("Q_def moment at right: %f in*lb\n", Qdef(4))

xvec1 = linspace(0,L/2);
xvec2 = linspace(L/2,L);
figure
hold on
plot(xvec1, subs(w1, xbar, xvec1))
plot(xvec2, subs(w2, xbar, xvec2-L/2))
title("Problem 6 w(x)")
xlabel("x (in)")
ylabel("Vertical displacement (in)")
hold off
figure
hold on
plot(xvec1, subs(sigma1, xbar, xvec1))
plot(xvec2, subs(sigma2, xbar, xvec2-L/2))
title("Problem 6 \sigma(x)")
xlabel("x (in)")
ylabel("Bending Stress (psi)")
hold off

fprintf("vertical dispalacement for 0<=x<=L/2 is w(x)=%s\n", subs(w1,xbar,x))
fprintf("vertical dispalacement for L/2<=x<=L is w(x)=%s\n", simplify(subs(w2,xbar,x-L/2)))
fprintf("bending stress for 0<=x<=L/2 is sigma(x)=%s\n", subs(sigma1,xbar,x))
fprintf("bending stress for L/2<=x<=L is sigma(x)=%s\n", simplify(subs(sigma2,xbar,x-L/2)))
