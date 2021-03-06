clear; close all; clc;
load('given.mat')
A = Ebar\Abar
B = Ebar\Bbar
C =  [0 0 0 1;
      0 1 0 0]
% C = eye(4)
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
fprintf("\tTherefore is controllable\n\n")

% Observability with only pitch rate feedback
Cq = [0 1 0 0];
fprintf("C for pitch rate only is \n")
disp(Cq)
O1 = Cq;
O2 = Cq*A;
O3 = Cq*A*A;
O4 = Cq*A*A*A;
O = [O1; O2; O3; O4];
fprintf("Observability with only pitch rate feedback\n")
disp(O)
fprintf("\tRank of observability should be %i and is %i\n", [rank(A), rank(O)])
fprintf("\tTherefore is observable\n\n")

% Observability with only pitch feedback
Ctheta = [0 0 0 1];
fprintf("C for pitch only is \n")
disp(Ctheta)
O1 = Ctheta;
O2 = Ctheta*A;
O3 = Ctheta*A*A;
O4 = Ctheta*A*A*A;
O = [O1; O2; O3; O4];
fprintf("Observability with only pitch feedback\n")
disp(O)
fprintf("\tRank of observability should be %i and is %i\n", [rank(A), rank(O)])
fprintf("\tTherefore is observable\n\n")

sso = ss(A, B, C, 0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
poles = pole(sso);
figure()
pzplot(sso)
title("Pole-Zero Map Open Loop")
fprintf("Open Loop System Stability:\n\tOpen Loop Poles\n")
disp(poles)
fprintf("\tAll real parts of poles are negative so open loop system is stable\n\n")
figure()
opt = stepDataOptions;
opt.StepAmplitude = deg2rad(10) ;
step(sso,opt)
title("Open Loop 10 Degree Step Response")
figure()
initial(sso, [0, 0, 0, deg2rad(10)])
title("Open Loop Response to 10 Degree Initial Pitch")

%% Part B
fprintf("***************Part B***************\n")
wn_sp = 16;
dr_sp = .8;
n = -wn_sp*(dr_sp);
w = wn_sp*sqrt(1-(dr_sp)^2);
lam1 = n + 1i*w;
lam2 = conj(lam1);
lambda = [lam1, lam2, poles(3), poles(4)];
Kack = acker(A, B, lambda);
% Kack = Ackermann(A, B, lambda);
fprintf("Ackermann's Formula\n")
fprintf("Desired Short Period Poles:\n")
disp([lam1; lam2])
fprintf("Gains are K = [%f  %f  %f  %f] \n", Kack)
fprintf("Poles using Ackermann gains:\n")
disp(eig(A-B*Kack))
ssack = ss(A-B*Kack, zeros(4,1), C, 0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
figure()
initial(ssack, [0, 0, 0, deg2rad(10)])
title("Ackermann gains Response to 10 Degree Initial Pitch")

[ol_vec, ol_lam] = eig(A);
% new_vec1 = [1+i*-1; -1+i*1; 0; 0];
% new_vec3 = [0; 0; 1+i*-1; -1+i*1];
% new_vec2 = conj(new_vec1);
% new_vec4 = conj(new_vec2);
% new_vec = [new_vec1, new_vec2, new_vec3, new_vec4];
% new_vec = [1-i, 1+i; 
%            -1+i, -1-i];
% new_vec = ol_vec([1,2],1:2)
% O = zeros(2,1); 
% D = [1 0 0 0;
%      0 1 0 0;];

new_vec = ol_vec(:,1:2);
O = zeros(4,1);
D = eye(4);
I4 = eye(4);
M1 = [(lambda(1)*I4-A), B; D, O];
M1inv = pinv(M1);
M2 = [(lambda(2)*I4-A), B; D, O];
M2inv = pinv(M2);
top_zeros = zeros(4,2);
des = [top_zeros; new_vec];
final1 = M1inv*des(:,1);
v1 = final1(1:4);
u1 = final1(5);
final2 = M2inv*des(:,2);
v2 = final2(1:4);
u2 = final2(5);
U = [u1, u2];
V = [v1, v2];
Kstruc = U*inv(C*V);

% Kstruc = EigStructAssign(A, B, C, lambda(1:2), new_vec);
ssc = ss(A-B*Kstruc*C,B,C,0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
figure()
initial(ssc, [0, 0, 0, deg2rad(10)])
title("Eigenstructure Assignment gains Response to 10 Degree Initial Pitch")

vd1 = [eye(4), zeros(4,1)];
fprintf("Eigenstructure Assignment using open loop eigen vectors\n")
fprintf("Desired Short Period Poles:\n")
disp([lam1; lam2])
fprintf("Gains are K = [%f  %f] \n", Kstruc)
fprintf("Poles using Eigenstructure Assignment gains:\n")
[es_vec, es_lam] = eig(A-B*Kstruc*C);
disp(es_lam)
fprintf("Achievable Eigen Vectors using Eigenstructure Assignment gains:\n")
disp(es_vec)


new_vec = [ol_vec(1,1), ol_vec(1,2)];
O = zeros(1,1); 
D = [1 0 0 0];
I4 = eye(4);
M1 = [(lambda(1)*I4-A), B; D, O];
M1inv = inv(M1);
M2 = [(lambda(2)*I4-A), B; D, O];
M2inv = inv(M2);
top_zeros = zeros(4,2);
des = [top_zeros; new_vec];
final1 = M1inv*des(:,1);
v1 = final1(1:4);
u1 = final1(5);
final2 = M2inv*des(:,2);
v2 = final2(1:4);
u2 = final2(5);
U = [u1, u2];
V = [v1, v2];
Kstruc = U*inv(C*V);

% Kstruc = EigStructAssign(A, B, C, lambda(1:2), new_vec);
ssc = ss(A-B*Kstruc*C,B,C,0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
figure()
initial(ssc, [0, 0, 0, deg2rad(10)])
title("Eigenstructure Assignment gains Response to 10 Degree Initial Pitch")

vd1 = [eye(4), zeros(4,1)];
fprintf("Eigenstructure Assignment of only alpha \n")
fprintf("Desired Short Period Poles:\n")
disp([lam1; lam2])
fprintf("Gains are K = [%f  %f] \n", Kstruc)
fprintf("Poles using Eigenstructure Assignment gains:\n")
[es_vec, es_lam] = eig(A-B*Kstruc*C);
disp(es_lam)
fprintf("Achievable Eigen Vectors using Eigenstructure Assignment gains:\n")
disp(es_vec)
fprintf("I am trying to change it as little as possible but for some reason\n" + ...
    "a pole I am not trying to change is going unstable. With just two \n" + ...
    "outputs only two sets of eigenvectors and values should be set\n")

%% Part C
fprintf("***************Part C***************\n")
alpha_max = 5;
q_max = 10;
VT_max = 100;
theta_max = 10;
demax = 10;
Q = diag([1/alpha_max^2, 1/q_max^2, 1/VT_max^2, 1/theta_max^2]);
R = 1/demax^2;
% R = 1e-16;
[Klqr, Slqr, CLPlqr] = lqr(sso, Q, R);
ssc = ss(A-B*Klqr,B,C,0, ...
    'StateName', {'alpha', 'q', 'VT', 'theta'}, ...
    'InputName', {'elevator angle'}, ...
    'OutputName', {'pitch', 'pitch rate'});
disp("Solution to Riccati")
disp(Slqr)
disp("Gains")
disp(Klqr)
figure()
initial(ssc, [0, 0, 0, deg2rad(10)])
title("LQR Response to 10 Degree Initial Pitch")

inp1deg = deg2rad(randn(201,1));
figure()
lsim(ssc, inp1deg, linspace(0, 10, 201), [0 0 0 deg2rad(10)])
title("LQR Response to 10 Degree Initial Pitch with noise")


%% Part D
fprintf("***************Part D***************\n")
Umat = [B, eye(4)];
sys2 = ss(A, Umat, C, 0);

inp_1deg = deg2rad(randn(201,1)*.1);
w_1 = [];
for ii = 1:length(inp_1deg)
    w_1(:,ii) = B.*inp_1deg(ii);
end
Qn_1 = cov(w_1');
% Qn = eye(4);
Rn_1 = cov(inp_1deg)*.1;
[Kest_1, L, P] = kalman(sys2, Qn_1, Rn_1);
reg_1 = lqgreg(Kest_1, Klqr);
closed_loop_1 = feedback(sso, -reg_1);
disp("Solution with noise on input of noise with normal distribution \sigma=.1 degree")
disp("standard deviation of .1 degree")
disp("Solution to Observer Riccati")
disp(P)
disp("Observer Gains")
disp(L)
figure()
initial(closed_loop, [0 0 0 deg2rad(10) 0 0 0 0], linspace(0, 30, 301))
title("LQG Response to 10 Degree Initial Pitch without noise")
figure()
lsim(closed_loop_1, inp_1deg, linspace(0, 30, 201), [0 0 0 deg2rad(10) 0 0 0 0])
title("LQG Response to 10 Degree Initial Pitch with input noise \sigma=.1 degree")


w = [];
for ii = 1:length(inp1deg)
    w(:,ii) = B.*inp1deg(ii);
end
Qn = cov(w');
% Qn = eye(4);
Rn = cov(inp1deg)*.1;
[Kest, L, P] = kalman(sys2, Qn, Rn);
reg = lqgreg(Kest, Klqr);
% disp(eig(reg))
closed_loop = feedback(sso, -reg);

disp("Solution with noise on input of noise with normal distribution \sigma=1 degree")
disp("standard deviation of 1 degree")
disp("Solution to Observer Riccati")
disp(P)
disp("Observer Gains")
disp(L)

figure()
lsim(closed_loop, inp1deg, linspace(0, 30, 201), [0 0 0 deg2rad(10) 0 0 0 0])
title("LQG Response to 10 Degree Initial Pitch with input noise \sigma=1 degree")


inp10deg = deg2rad(randn(201,1)*10);
w10 = [];
for ii = 1:length(inp10deg)
    w10(:,ii) = B.*inp10deg(ii);
end
Qn10 = cov(w10');
% Qn = eye(4);
Rn10 = cov(inp10deg);
[Kest10, L, P] = kalman(sys2, Qn10, Rn10);
reg10 = lqgreg(Kest10, Klqr);
closed_loop10 = feedback(sso, -reg10);
disp("Solution with noise on input of noise with normal distribution \sigma=10 degree")
disp("standard deviation of 10 degree")
disp("Solution to Observer Riccati")
disp(P)
disp("Observer Gains")
disp(L)
figure()
lsim(closed_loop10, inp10deg, linspace(0, 30, 201), [0 0 0 deg2rad(10) 0 0 0 0])
title("LQG Response to 10 Degree Initial Pitch with input noise \sigma=10 degree")

disp("After testing with different magnitudes of noise it is clear that the")
disp("the obvious conclusion is correct that the systems work better with ")
disp("less noise. It seems like 10 degree standard deviation for the noise")
disp("is far too much. The pitch rate seems more effected than the pitch by ")
disp("the noise for some reason. Noise was modeled as a noisy input in lsim ")
disp("and it was scaled to 10, 1, and .1 degrees. The input was shifted to  ")
disp("radians before being input to lsim")

% %% Functions
% 
% function K = Ackermann(A, B, lambda)
%     n = length(A);
% 
%     syms s
%     funs = 1;
%     for ii = 1:length(lambda)
%         funs = funs*(s-lambda(ii));
%     end
%     fdes_s = collect(funs);
%     coef = sym2poly(fdes_s);
%     fdes_A = 0;
%     m = length(coef);
%     for ii = 1:m
%         fdes_A = fdes_A + coef(ii)*A^(n-1);
%     end
% %     fdes_A = coef(1).*A*A + coef(2).*A + coef(3).*eye(n);
% 
%     Con = [];
%     for ii = 1:n
%         Con = [Con, A^(ii-1)*B];
%     end
% 
%     leftv = [zeros(1,n-1), 1];
%     K = leftv*inv(Con)*fdes_A;
% %     fprintf("Control Gains are K = \n")
% %     disp(K)
% end
% 
% function K = EigStructAssign(A, B, C, des_lambda, des_vec)
%     n = length(A);
%     m = size(B,2);
%     p = size(des_vec,1);
% 
%     I = eye(n);
%     D = diag(ones(1,p));
%     O = zeros(p,(n+m)-p);
%     top = zeros(n,1);
%     uv = @(lambda,vec) pinv([lambda.*I - A, B; D, O])*[top;vec];
% 
%     for ii = 1:length(des_lambda)
%         uvii = uv(des_lambda(ii), des_vec(:,ii));
%         vd(:,ii) = uvii(1:n);
%         ud(:,ii) = uvii(n+1:end);
%     end
%     K = ud*pinv(C*vd);
%     fprintf("Control Gains are K = \n")
%     disp(K)
%     poles = eig(A-B*K*C);
%     fprintf("Closed Loop Poles are \\lamda = \n")
%     disp(poles)
%     fprintf("Desired Eigen Vectors are v = \n")
%     disp(des_vec)
%     fprintf("Achievable Eigen Vectors are v = \n")
%     disp(vd)
% end
% 
% function K = outputLQR(A, B, C, Q, R, x0, tol)
%     n = length(A);
%     m = size(B,2);
%     p = size(C,1);
%     P = sym("P", [n n]);
%     S = sym("S", [n n]);
% %     syms P S
%     if x0 == 0
%         x0 = zeros(n,1);
%     end
%     X = x0*x0';
%     X = eye(n);
%     kk = 0;
%     Kk = zeros(m,p);
%     Ak = A - B*Kk*C;
%     eqPk = 0 == Ak.'*P + P*Ak + C.'*Kk.'*R*Kk*C + Q;
%     eqSk = 0 == Ak*S + S*S.' + X;
%     sol = vpasolve([eqSk],[S]);
% %     Pk = sol.Pk
%     Sk = zeros(n);
%     Pk = zeros(n);
%     Sk1 = zeros(n);
%     Pk1 = zeros(n);
%     for ii = 1:n
%         for jj = 1:n
%             Sk(ii,jj) = eval(['sol.S', int2str(ii), '_', int2str(jj)]);
%             Pk(ii,jj) = eval(['sol.P', int2str(ii), '_', int2str(jj)]);
%         end
%     end
%     Jk = (1/2)*trace(Pk*X);
%     alpha = 1e-2;
%     err  = 1;
%     while err > tol
%         Ak = A - B*Kk*C;
%         eqPk = 0 == Ak'*P + P*Ak + C'*Kk'*R*Kk*C + Q;
%         eqSk = 0 == Ak*S + S*S' + X;
%         sol = vpasolve([eqPk,eqSk],[S,P]);
%         for ii = 1:n
%             for jj = 1:n
%                 Sk1(ii,jj) = eval(['sol.S', int2str(ii), '_', int2str(jj)]);
%                 Pk1(ii,jj) = eval(['sol.P', int2str(ii), '_', int2str(jj)]);
%             end
%         end
%         Jk1 = (1/2)*trace(Pk1*X);
%         deltaK = inv(R)*B'*Pk1*Sk1*C'*inv(C*Sk1*C')-Kk;
%         Kk1 = Kk + alpha*deltaK;
%         err = abs(Jk1-Jk);
%         kk = kk + 1;
%         Sk = Sk1;
%         Pk = Pk1;
%         Jk = Jk1;
%         Kk = Kk1;
%     end
%     K = Kk;
% end