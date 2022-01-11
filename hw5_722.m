%% Problem 1
clear; close all; clc;
syms x;

Ks = 5/6;
E = 30*10^6;
b = 3;
c = 2;
L = 144;
h = L/2;
I = b*(c*2)^3/12;
EI = E*I;
G = 11*10^6;
A = b*c*2;
GAKs = G*A*Ks;

ksp = 2000;
q0 = 50;
M0 = 10000;
F0 = 2000;

quad1 = (1-x/h)*(1-2*x/h);
quad2 = (4*x/h)*(1-x/h);
quad3 = (-x/h)*(1-2*x/h);

lin1 = 1-x/h;
lin2 = x/h;

LAMBDA = EI/(GAKs*h);
mu0 = 12*LAMBDA;

Ke = (2*EI/(mu0*h^3))*[6, -3*h, -6, -3*h;
    -3*h, h^2*(1.5+6*LAMBDA), 3*h, h^2*(1.5-6*LAMBDA);
    -6, 3*h, 6, 3*h;
    -3*h, h^2*(1.5-6*LAMBDA), 3*h, h^2*(1.5+6*LAMBDA)];
K1 = Ke;
K2 = Ke;
Kg = zeros(6);
for ii = 1:4
    for jj = 1:4
        Kg(ii,jj) = Kg(ii,jj) + K1(ii,jj);
        Kg(ii+2,jj+2) = Kg(ii+2,jj+2) + K1(ii,jj);
    end
end
Kg(5,5) = Kg(5,5) + ksp;
q = q0*h/12*[6;-h;6;h;0;0];
Qb = [M0; F0; 0];
DELTAt = [0;0;0];
DELTAb = inv(Kg(4:6,4:6))*(q(4:6)+Qb-Kg(4:6,1:3)*DELTAt);
Qt = Kg(1:3,1:3)*DELTAt + Kg(1:3,4:6)*DELTAb;
DELTA = [DELTAt; DELTAb];
qc = eval(int(quad2*q0,0,h));
wc1 = (6*h/(32*GAKs))*(qc)+((DELTA(1)+DELTA(3))/2)+h*((DELTA(4)-DELTA(2))/2);
wc2 = (6*h/(32*GAKs))*(0)+((DELTA(3)+DELTA(5))/2)+h*((DELTA(6)-DELTA(4))/2);
w1 = quad1*DELTA(1) + quad2*wc1 + quad3*DELTA(3);
w2 = quad1*DELTA(3) + quad2*wc2 + quad3*DELTA(5);
psi1 = lin1*DELTA(2) + lin2*DELTA(4);
psi2 = lin1*DELTA(4) + lin2*DELTA(6);

V1 = GAKs*(psi1+diff(w1));
V2 = GAKs*(psi1+diff(w2));
M1 = EI*diff(psi1);
M2 = EI*diff(psi2);
sigma1 = E*c*diff(psi1);
sigma2 = E*c*diff(psi2);

xvec1 = linspace(0,h,1e2);
xvec2 = linspace(h,L,1e2);

figure
plot(xvec1, subs(w1,x,xvec1), 'b')
plot(xvec2, subs(w2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Vertical Displacement (in)")
title("Problem 1")

figure
plot(xvec1, subs(sigma1,x,xvec1), 'b')
plot(xvec2, subs(sigma2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Bending Stress (psi)")
title("Problem 1")

fprintf("Problem 1 \n")
fprintf("Q1_def = %f lb\n",-eval(subs(V1,x,0)))
fprintf("Q1_equ = %f lb\n",Qt(1))
fprintf("Q2_def = %f in-lb\n",-eval(subs(M1,x,0)))
fprintf("Q2_equ = %f in-lb\n",Qt(2))
fprintf("Q5_def = %f lb\n",eval(subs(V2,x,h)))
fprintf("Q5_equ = %f lb\n",Qb(2))
fprintf("Q6_def = %f in-lb\n",eval(subs(M2,x,h)))
fprintf("Q6_equ = %f in-lb\n",Qb(3))

%% Problem 3
fprintf("Problem 3 \n")
clear;
syms x h;

Ks = 5/6;
E = 30*10^6;
b = 3;
c = 2;
L1 = 50*sqrt(2);
L2 = 72;
I = b*(c*2)^3/12;
EI = E*I;
G = 11*10^6;
A = b*c*2;
GAKs = G*A*Ks;

q0 = 50;

alpha = -45*(pi/180);

% cubic1 = 1-3*(x/h)^2+2*(x/h)^3;
% cubic2 = -x*(1-x/h)^2;
% cubic3 = 3*(x/h)^2-2*(x/h)^3;
% cubic4 = -x*((1-x/h)^2-(x/h));

quad1 = (1-x/h)*(1-2*x/h);
quad2 = (4*x/h)*(1-x/h);
quad3 = (-x/h)*(1-2*x/h);

lin1 = 1-x/h;
lin2 = x/h;

LAMBDA = EI/(GAKs*h);
mu0 = 12*LAMBDA;
mu = A*mu0*h^2/(2*I);

Ke = (2*EI/(mu0*h^3))*[mu, 0, 0,-mu, 0, 0;
    0, 6, -3*h, 0, -6, -3*h;
    0, -3*h, h^2*(1.5+6*LAMBDA), 0, 3*h, h^2*(1.5-6*LAMBDA);
    -mu, 0, 0, mu, 0, 0;
    0, -6, 3*h, 0, 6, 3*h;
    0, -3*h, h^2*(1.5-6*LAMBDA), 0, 3*h, h^2*(1.5+6*LAMBDA)];
K1loc = subs(Ke,h,L1);
K2 = subs(Ke,h,L2);
T = [cos(alpha), sin(alpha), 0, 0, 0, 0;
    -sin(alpha), cos(alpha), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(alpha), sin(alpha), 0;
    0, 0, 0, -sin(alpha), cos(alpha), 0;
    0, 0, 0, 0, 0, 1];
K1 = T'*K1loc*T;
Kg = zeros(9);
for ii = 1:6
    for jj = 1:6
        Kg(ii,jj) = Kg(ii,jj) + K1(ii,jj);
        Kg(ii+2,jj+2) = Kg(ii+3,jj+3) + K1(ii,jj);
    end
end
q = q0*L2/12*[0;0;0;0;6;-L2;0;6;L2];
Qb = [0; 0; 0];
DELTAa = [0;0;0];
DELTAc = DELTAa;
DELTAb = inv(Kg(4:6,4:6))*(q(4:6)+Qb-Kg(4:6,1:3)*DELTAa-Kg(4:6,7:9)*DELTAc);
Qa = Kg(1:3,1:3)*DELTAa + Kg(1:3,4:6)*DELTAb + Kg(1:3,7:9)*DELTAc - q(1:3);
Qc = Kg(7:9,1:3)*DELTAa + Kg(7:9,4:6)*DELTAb + Kg(7:9,7:9)*DELTAc - q(7:9);
DELTA = [DELTAa; DELTAb; DELTAc];

DELTA1_loc = T*DELTA(1:6);
qc = eval(int(subs(quad2,h,L2)*q0,0,L2));
wc1 = (6*L1/(32*GAKs))*(0)+((DELTA(2)+DELTA(5))/2)+L1*((DELTA(6)-DELTA(3))/2);
wc2 = (6*L2/(32*GAKs))*(qc)+((DELTA(5)+DELTA(8))/2)+L2*((DELTA(9)-DELTA(6))/2);
w1 = subs(quad1,h,L1)*DELTA(2) + subs(quad2,h,L1)*wc1 + subs(quad3,h,L1)*DELTA(5);
w2 = subs(quad1,h,L2)*DELTA(5) + subs(quad2,h,L2)*wc2 + subs(quad3,h,L2)*DELTA(8);
psi1 = subs(lin1,h,L1)*DELTA(3) + subs(lin2,h,L1)*DELTA(6);
psi2 = subs(lin1,h,L2)*DELTA(6) + subs(lin2,h,L2)*DELTA(9);
u1 = subs(lin1,h,L1)*DELTA(1) + subs(lin2,h,L1)*DELTA(4);
u2 = subs(lin1,h,L2)*DELTA(4) + subs(lin2,h,L2)*DELTA(7);

V1 = GAKs*(psi1+diff(w1));
V2 = GAKs*(psi1+diff(w2));
M1 = EI*diff(psi1);
M2 = EI*diff(psi2);
sigma1 = -E*c*diff(psi1)+E*diff(u1);
sigma2 = -E*c*diff(psi2)+E*diff(u2);

xvec1 = linspace(0,L1,1e2);
xvec2 = linspace(L1,L2,1e2);

figure
plot(xvec1, subs(w1,x,xvec1), 'b')
plot(xvec2, subs(w2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Transverse Displacement (in)")
title("Problem 3: Timoshenko")

figure
plot(xvec1, subs(sigma1,x,xvec1), 'b')
plot(xvec2, subs(sigma2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Max Axial Stress (psi)")
title("Problem 3: Timoshenko")

%% Problem 2
fprintf("Problem 2 \n")
clear;
syms x h;

Ks = 5/6;
E = 30*10^6;
b = 3;
c = 2;
L1 = 50*sqrt(2);
L2 = 72;
I = b*(c*2)^3/12;
EI = E*I;
G = 11*10^6;
A = b*c*2;
GAKs = G*A*Ks;

q0 = 50;

alpha = -45*(pi/180);

cubic1 = 1-3*(x/h)^2+2*(x/h)^3;
cubic2 = -x*(1-x/h)^2;
cubic3 = 3*(x/h)^2-2*(x/h)^3;
cubic4 = -x*((1-x/h)^2-(x/h));

% quad1 = (1-x/h)*(1-2*x/h);
% quad2 = (4*x/h)*(1-x/h);
% quad3 = (-x/h)*(1-2*x/h);

lin1 = 1-x/h;
lin2 = x/h;

LAMBDA = EI/(GAKs*h);
mu0 = 12*LAMBDA;
mu = A*mu0*h^2/(2*I);

Ke = (2*EI/(h^3))*[mu, 0, 0,-mu, 0, 0;
    0, 6, -3*h, 0, -6, -3*h;
    0, -3*h, 2*h^2, 0, 3*h, h^2;
    -mu, 0, 0, mu, 0, 0;
    0, -6, 3*h, 0, 6, 3*h;
    0, -3*h, h^2, 0, 3*h, 2*h^2];
K1loc = subs(Ke,h,L1);
K2 = subs(Ke,h,L2);
T = [cos(alpha), sin(alpha), 0, 0, 0, 0;
    -sin(alpha), cos(alpha), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(alpha), sin(alpha), 0;
    0, 0, 0, -sin(alpha), cos(alpha), 0;
    0, 0, 0, 0, 0, 1];
K1 = T'*K1loc*T;
Kg = zeros(9);
for ii = 1:6
    for jj = 1:6
        Kg(ii,jj) = Kg(ii,jj) + K1(ii,jj);
        Kg(ii+2,jj+2) = Kg(ii+3,jj+3) + K1(ii,jj);
    end
end
q = q0*L2/12*[0;0;0;0;6;-L2;0;6;L2];
Qb = [0; 0; 0];
DELTAa = [0;0;0];
DELTAc = DELTAa;
DELTAb = inv(Kg(4:6,4:6))*(q(4:6)+Qb-Kg(4:6,1:3)*DELTAa-Kg(4:6,7:9)*DELTAc);
Qa = Kg(1:3,1:3)*DELTAa + Kg(1:3,4:6)*DELTAb + Kg(1:3,7:9)*DELTAc - q(1:3);
Qc = Kg(7:9,1:3)*DELTAa + Kg(7:9,4:6)*DELTAb + Kg(7:9,7:9)*DELTAc - q(7:9);
DELTA = [DELTAa; DELTAb; DELTAc];

DELTA1_loc = T*DELTA(1:6);
w1 = subs(cubic1,h,L1)*DELTA(2) + subs(cubic2,h,L1)*DELTA(3) + subs(cubic3,h,L1)*DELTA(5) + subs(cubic4,h,L1)*DELTA(6);
w2 = subs(cubic1,h,L2)*DELTA(5) + subs(cubic2,h,L2)*DELTA(6) + subs(cubic3,h,L2)*DELTA(8) + subs(cubic4,h,L2)*DELTA(9);
u1 = subs(lin1,h,L1)*DELTA(1) + subs(lin2,h,L1)*DELTA(4);
u2 = subs(lin1,h,L2)*DELTA(4) + subs(lin2,h,L2)*DELTA(7);

sigma1 = -EI*c*diff(diff(w1))+E*diff(u1);
sigma2 = -EI*c*diff(diff(w2))+E*diff(u2);

xvec1 = linspace(0,L1,1e2);
xvec2 = linspace(L1,L2,1e2);

figure
plot(xvec1, subs(w1,x,xvec1), 'b')
plot(xvec2, subs(w2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Transverse Displacement (in)")
title("Problem 2: Euler Bernoulli")

figure
plot(xvec1, subs(sigma1,x,xvec1), 'b')
plot(xvec2, subs(sigma2,x,xvec1), 'b')
xlabel("x (in)")
ylabel("Max Axial Stress (psi)")
title("Problem 2: Euler Bernoulli")

%% Problem 4
fprintf("Problem 4 \n")
clear;
syms x theta2 theta3 w3 lambda;

E = 30*10^6;
b = 3;
c = 2;
L = 144;
h = L/2;
I = b*(c*2)^3/12;
EI = E*I;
rho = 9*10^-3;
A = b*2*c;

ksp = 2000;

cubic1 = 1-3*(x/h)^2+2*(x/h)^3;
cubic2 = -x*(1-x/h)^2;
cubic3 = 3*(x/h)^2-2*(x/h)^3;
cubic4 = -x*((1-x/h)^2-(x/h));

Ke = (2*EI/(h^3))*[6, -3*h, -6, -3*h;
    -3*h, 2*h^2, 3*h, h^2;
    -6, 3*h, 6, 3*h;
    -3*h, h^2, 3*h, 2*h^2];
Me = (rho*A*h/420)*[156, -22*h, 54, 13*h;
                    -22*h, 4*h^2, -13*h, -3*h^2;
                    54, -13*h, 156, 22*h;
                    13*h, -3*h, 22*h, 4*h] + ...
     (rho*I/(30*h))*[36, -3*h, -36, -3*h;
                     -3*h, 4*h^2, 3*h, -h^2;
                     -36, 3*h, 36, 3*h;
                     -3*h, -h^2, 3*h, 4*h^2];
K1 = Ke;
K2 = Ke;
M1 = Me;
M2 = Me;
Kg = zeros(6);
Mg = zeros(6);
for ii = 1:4
    for jj = 1:4
        Kg(ii,jj) = Kg(ii,jj) + K1(ii,jj);
        Kg(ii+2,jj+2) = Kg(ii+2,jj+2) + K2(ii,jj);
        Mg(ii,jj) = Mg(ii,jj) + M1(ii,jj);
        Mg(ii+2,jj+2) = Mg(ii+2,jj+2) + M2(ii,jj);
    end
end
Kg(5,5) = Kg(5,5) + ksp;
eigenvalues = vpasolve(0==det(Kg(4:6,4:6)-lambda*Mg(4:6,4:6)),lambda);
fprintf("Fundamental frequencies are: %f, %f, and %f%+fj \n",sqrt(eigenvalues(2)), sqrt(eigenvalues(3)), real(sqrt(eigenvalues(1))), imag(sqrt(eigenvalues(1))))


