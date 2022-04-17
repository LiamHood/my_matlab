clear; close all; clc;
load('given.mat')
A = Ebar\Abar
B = Ebar\Bbar
C =  [0 1 0 0; 0 0 0 1]
% x = [alpha, q, VT, theta]
% angle of attack, pitch rate, true air speed, pitch angle

%% Part A
fprintf("***************Part A***************\n")
% Controllability
U1 = B;
U2 = A*B;
U3 = A*A*B;
U4 = A*A*A*B;
U = [U1 U2 U3 U4];
fprintf("Controllability\n")
disp(U)
fprintf("\tRank of controllability should be %i and is %i\n", [rank(A), rank(U)])
fprintf("\tTherefore is controllable\n")

% Observability with only pitch rate feedback
Cq = [0 1 0 0];
O1 = Cq;
O2 = Cq*A;
O3 = Cq*A*A;
O4 = Cq*A*A*A;
O = [O1; O2; O3; O4];
fprintf("Observability with only pitch rate feedback\n")
disp(O)
fprintf("\tRank of observability should be %i and is %i\n", [rank(A), rank(O)])
fprintf("\tTherefore is observable\n")

% Observability with only pitch feedback
Ctheta = [0 0 0 1];
O1 = Ctheta;
O2 = Ctheta*A;
O3 = Ctheta*A*A;
O4 = Ctheta*A*A*A;
O = [O1; O2; O3; O4];
fprintf("Observability with only pitch feedback\n")
disp(O)
fprintf("\tRank of observability should be %i and is %i\n", [rank(A), rank(O)])
fprintf("\tTherefore is observable\n")

sso = ss(A, B, C, 0);
poles = pole(sso);
figure()
pzplot(sso)
fprintf("Open Loop System Stability:\n\tOpen Loop Poles\n")
disp(poles)
fprintf("\tAll real parts of poles are negative so open loop system is stable\n")
figure()
opt = stepDataOptions;
opt.StepAmplitude = -1.2;
step(sso,opt)
figure()
initial(sso, [0, 0, 0, 10])

%% Part B
fprintf("***************Part B***************\n")
wn_sp = 16;
dr_sp = .8;
lam1 = -dr_sp + 1i*wn_sp;
lam2 = conj(lam1);
lambda = [lam1, lam2, poles(3), poles(4)];
Kack = Ackermann(A, B, lambda);
fprintf("Ackermann's Formula\n")
fprintf("\t Gains are K = [%f  %f  %f  %f] \n", Kack)

[ol_vec,ol_lam] = eig(A);
new_vec1 = [1+i*-1; -1+i*1; 0; 0]
new_vec3 = [0; 0; 1+i*-1; -1+i*1]
new_vec2 = conj(new_vec1)
new_vec4 = conj(new_vec2)
new_vec = [new_vec1, new_vec2, new_vec3, new_vec4]
Kstruc = EigStructAssign(A, B, C, lambda, new_vec);

%% Part C
fprintf("***************Part C***************\n")
alpha_max = 5;
q_max = 10;
VT_max = 100;
theta_max = 20;
demax = 20;
Q = diag([1/alpha_max^2, 1/q_max^2, 0, 1/theta_max^2]);
R = 1/demax^2;
[Klqr, Slqr, CLPlqr] = lqr(sso, Q, R);
ssc = ss(A-B*Klqr,zeros(length(A),1),C,0);
disp("Solution to Riccati")
disp(Slqr)
disp("Gains")
disp(Klqr)
figure(1)
initial(ssc, [0, 0, 0, 10])

%% Functions

function K = Ackermann(A, B, lambda)
    n = length(A);

    syms s
    funs = 1;
    for ii = 1:length(lambda)
        funs = funs*(s-lambda(ii));
    end
    fdes_s = collect(funs);
    coef = sym2poly(fdes_s);
    fdes_A = 0;
    m = length(coef);
    for ii = 1:m
        fdes_A = fdes_A + coef(ii)*A^(n-1);
    end
%     fdes_A = coef(1).*A*A + coef(2).*A + coef(3).*eye(n);

    Con = [];
    for ii = 1:n
        Con = [Con, A^(ii-1)*B];
    end

    leftv = [zeros(1,n-1), 1];
    K = leftv*inv(Con)*fdes_A;
%     fprintf("Control Gains are K = \n")
%     disp(K)
end

function K = EigStructAssign(A, B, C, des_lambda, des_vec)
    poles = eig(A);
%     fprintf("Open Loop Poles are \\lamda = \n")
    disp(poles)
    
    n = length(A);

    Con = [];
    for ii = 1:n
        Con = [Con, A^(ii-1)*B];
    end
%     fprintf("Controllable since controllability matrix is rank %i and A is rank %i\n\n", [rank(Con), rank(A)])

    
    m = size(B,2);
    p = length(des_lambda);

    I = eye(n);
    D = diag(ones(1,p));
    O = zeros(p,(n+m)-p);
    top = zeros(n,1);
    uv = @(lambda,vec) pinv([lambda.*I - A, B; D, O])*[top;vec];

    for ii = 1:p
        uvii = uv(des_lambda(ii), des_vec(:,ii));
        vd(:,ii) = uvii(1:n);
        ud(:,ii) = uvii(n+1:end);
    end
    K = ud*pinv(C*vd);
    fprintf("Control Gains are K = \n")
    disp(K)
    poles = eig(A-B*K*C);
    fprintf("Closed Loop Poles are \\lamda = \n")
    disp(poles)
    fprintf("Desired Eigen Vectors are v = \n")
    disp(des_vec)
    fprintf("Achievable Eigen Vectors are v = \n")
    disp(vd)
end
