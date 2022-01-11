function [M] = simGyroMotion(I_t, I_a, w_0, M_0, tspan)

H_i = [0; 0; 1];
I = diag([I_t, I_t, I_a]);
h_0 = I*w_0;

% Define initial PRECESSION angle
theta_0 = 0;
% Define initial NUTATION angle
%h_0 = I*
gamma_0 = acos( h_0(3,1)/norm(h_0) );
% Define initial SPIN angle
mu_0 = atan2(w_0(1), w_0(2));

beta = atan(I_a/I_t*tan(gamma_0));
alpha = abs(gamma_0 - beta);

% Draw Body Graphical Components
m = 1500;
r = sqrt(2*I_a/m);
l = sqrt((I_t - I_a/2)*12/m);

spaceCone = cone(cos(alpha), sin(alpha) ,sin(alpha), 60, [0 0 0], H_i(:,1), 'facealpha', .4);%psi_dot);
bodyCone = cone(cos(beta), sin(beta), sin(beta), 60, [0 0 0], [0; 0; 1], 'faceColor', 'red', 'facealpha', .4);

b1_vec = arrow([0 1], [0 0], [0 0], 'facecolor', 'red');
b2_vec = arrow([0 0], [0 1], [0 0], 'facecolor', 'red');
b3_vec = arrow([0 0], [0 0], [0 1], 'facecolor', 'red');
w_vec = arrow([0 w_0(1,1)/norm(w_0)], [0 w_0(2,1)/norm(w_0)], [0 w_0(3,1)/norm(w_0)], 'facecolor', 'yellow');
rb = cylin(r,r,l, 120, 'top', 'facecolor',[.1 .1 .1]);
w_vec_0 = arrow([0 w_0(1,1)/norm(w_0)], [0 w_0(2,1)/norm(w_0)], [0 w_0(3,1)/norm(w_0)], 'facecolor', 'green');
translate(rb, 0,0,-l/2);

% Draw inertial graphical componenets
i1_vec = arrow([0 1], [0 0], [0 0], 'facecolor', 'blue');
i2_vec = arrow([0 0], [0 1], [0 0], 'facecolor', 'blue');
i3_vec = arrow([0 0], [0 0], [0 1], 'facecolor', 'blue');

% Get everything into initial position
rotate([b1_vec, b2_vec, b3_vec rb w_vec w_vec_0 bodyCone], [0 0 1]', theta_0*180/pi, [0 0 0]');
rotate([b1_vec, b2_vec, b3_vec rb w_vec w_vec_0 bodyCone], Cz(theta_0)'*[1 0 0]', gamma_0*180/pi, [0 0 0]');
rotate([b1_vec, b2_vec, b3_vec rb w_vec w_vec_0 bodyCone], Cz(theta_0)'*Cx(gamma_0)'*[0 0 1]', mu_0*180/pi, [0 0 0]');

% Simulate the motion
display('Calculating...')

% Calculating more then we need...
[H, H_i, w, wi, b3, b3i, E, E_dot, q, t, E_0, E_dot_0] =...
    calc_dynamics_2(I, M_0, w_0, tspan);

M(1) = getframe;
% find the deltas
dtheta = diff(E(:,1));
dgamma = diff(E(:,2));
dmu = diff(E(:,3));

% loop through the results for the simulation
for i = 1:length(dtheta)
    rotate([bodyCone b1_vec, b2_vec, b3_vec rb w_vec],...
        [0 0 1]', dtheta(i)*180/pi, [0 0 0]');
    rotate([bodyCone b1_vec, b2_vec, b3_vec rb w_vec],...
        Cz(E(i,1))'*[1 0 0]', dgamma(i)*180/pi, [0 0 0]');
    rotate([bodyCone b1_vec, b2_vec, b3_vec rb w_vec],...
        Cz(E(i,1))'*Cx(E(i,2))'*[0 0 1]', dmu(i)*180/pi, [0 0 0]');
    M(i+1) = getframe;
    drawnow
end