%% A.3 Optimal Ascent
clear; close all; clc;

% define parameters
global g_accel Thrust2Weight    % pass to DE function
h = 185.2e3;    % meters, final altitude
vc = 1.627e3;   % m/s, circular speed at 100 nmi
g_accel = 1.62; % m/s^2, gravitational acceleration of the moon
Thrust2Weight = 3; % Thrust to weight of lunar ascent vehicle in lunar G's
rad2deg = 180/pi;

%% Boundary Conditions
global x0 y0 Vx0 Vy0 yf Vxf Vyf     % pass to BC function

% initial conditions
x0 = 0;         % meters
y0 = 0;
Vx0 = 0;        % meters/second
Vy0 = 0;

% final conditions
yf = h;         % meters, final altitude
Vxf = vc;       % m/s, final downrange velocity
Vyf = 0;        % final vertical velocity

%% Initial guess
t0 = 0;
% list initial conditions in yinit, use 0 if unknown
yinit = [x0 y0 Vx0 Vy0 0 0];        % guess for initial state and costate
tf_guess = 700;                     % sec, initial guess for final time

% in order to leave tf as a free variable to optimize, we parameterize by
% tf. tau = time/tf
Nt = 41; % number of time points
tau = linspace(0,1,Nt)';    % non-dimensional time vector

% use bvpinit to get the initial guess
% requires non-dimensional time vector; initial states, costates and 
% guesses; and guess for the final time
solinit = bvpinit(tau,yinit,tf_guess);

%% solution

% call bvp4c to solve TPBVP. Point solver to the functions containing the
% differential equations and the boundary conditions and provide it with
% the initial guess
sol = bvp4c(@ascent_odes_tf, @ascent_bcs_tf, solinit);

tf = sol.parameters(1);     % final time

% Evaluate at all times in the nondimensional time vector tau and store the
% state variable in matrix z
z = deval(sol,tau);

% Convert to dimensional time so that it means something
time = t0 + tau.*(tf-t0);

% extract sol for each state from matrix z
x_sol = z(1,:);
y_sol = z(2,:);
vx_sol = z(3,:);
vy_sol = z(4,:);
lambda2_bar_sol = z(5,:);
lambda4_bar_sol = z(6,:);

%% Plots

figure
subplot(3,2,1)
plot(time,x_sol/1000);
ylabel('x, km', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;

title('Optimal Ascent from Flat Moon')

subplot(3,2,2)
plot(time,y_sol/1000);
ylabel('y, km', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;

subplot(3,2,3)
plot(time,vx_sol/1000);
ylabel('V_x, km/s', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;

subplot(3,2,4)
plot(time,vy_sol/1000);
ylabel('V_y, km/s', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;

subplot(3,2,5)
plot(time,rad2deg*atan(lambda4_bar_sol));
ylabel('\alpha, deg', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;

subplot(3,2,6)
plot(time,lambda4_bar_sol);
ylabel('\lambda_4', 'fontsize', 14)
xlabel('Time, sec', 'fontsize', 14)
hold on; grid on;
hold off;

figure
plot(x_sol/1000, y_sol/1000)
grid on; axis equal;
xlabel('Downrange Position, x, km', 'fontsize', 14)
ylabel('Altitude, y, km', 'fontsize', 14)

%% functions

function PSI = ascent_bcs_tf(Y0, Yf, tf)
    global x0 y0 Vx0 Vy0 yf Vxf Vyf
    PSI = [ Y0(1) - x0
            Y0(2) - y0
            Y0(3) - Vx0
            Y0(4) - Vy0
            Yf(2) - yf
            Yf(3) - Vxf
            Yf(4) - Vyf
          ];
end

function dX_dtau = ascent_odes_tf(tau, X, tf)
    global g_accel Thrust2Weight
    Acc = Thrust2Weight*g_accel;    % Acceleration (F/m) of the Ascent vehicle m/s^2
    xdot = X(3);
    ydot = X(4);
    Vxdot = Acc*(1/sqrt(1+X(6)^2));
    Vydot = Acc*(X(6)/sqrt(1+X(6)^2)) - g_accel;
    lambda2_bar_dot = 0;
    lambda4_bar_dot = -X(5);
    dX_dtau = tf*[xdot; ydot; Vxdot; Vydot; lambda2_bar_dot; lambda4_bar_dot];
end