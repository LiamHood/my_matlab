function [H, Hi, w, wi, b3, b3i, E, E_dot, q, t, E_0, E_dot_0] = calc_dynamics_2(J, M_0, w_0, tspan)

H_0 = J*w_0; %initial angular momentum in body frame

% Initial PRECESSION angle
theta_0 = atan2(w_0(1), w_0(2));

% Initial NUTATION angle
gamma_0 = acos( H_0(3,1)/norm(H_0) );

% Initial SPIN angle
mu_0 = 0;

E_0 = [mu_0; gamma_0; theta_0];

E_dot_0 = LBI(mu_0, gamma_0, theta_0, 313)*w_0;

Q0 = Cz(mu_0)*Cx(gamma_0)*Cz(theta_0);
q_0(4,1) = sqrt(1+trace(Q0))/2;
q_0(1,1) = (Q0(3,2)-Q0(2,3))/4/q_0(4,1);
q_0(2,1) = (Q0(1,3)-Q0(3,1))/4/q_0(4,1);
q_0(3,1) = (Q0(2,1)-Q0(1,2))/4/q_0(4,1);


init_cond = [w_0; E_0; q_0];

opt = odeset('relTol', 1e-8, 'abstol', 1e-8);

[t, y] = ode45(@euler_fun, tspan, init_cond, opt, J, M_0);

% Body Rates
w = y(:,1:3);
% Precession
theta = y(:,4);
% Nutation
gamma = y(:,5);
% Spin
mu = y(:,6);

E = [theta, gamma, mu];

q = y(:,7:10);

b1 = [1; 0; 0];
b2 = [0; 1; 0];
b3 = [0; 0; 1];

wi = zeros(3, length(y));
H = zeros(3, length(y));
Hi = zeros(3, length(y));
b1i = zeros(3, length(y));
b2i = zeros(3, length(y));
b3i = zeros(3, length(y));
E_dot = zeros(3, length(y));

for i = 1:size(t,1)
    CIB = (Cz(mu(i))*Cx(gamma(i))*Cz(theta(i)))';
    wi(:,i) = CIB*w(i,:)';
    H(:,i) = J*w(i,:)';
    Hi(:,i) = CIB*H(:,i);
    b1i(:,i) = CIB*b1;
    b2i(:,i) = CIB*b2;
    b3i(:,i) = CIB*b3;
    E_dot(:,i) = LBI(theta(i), gamma(i), mu(i), 313)*w(i,:)';
end