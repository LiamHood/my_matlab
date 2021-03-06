%% B.4 Optimal Ascent
clear; close all; clc;

% define parameters
global g_accel Vc h eta beta f m0 mdot    % pass to DE function
h = 180e3;    % meters, final altitude
Vc = sqrt(3.9860044e5/(6378.14+h/1000))*1000;   % m/s, circular speed at 180 km
g_accel = 9.80665;  % m/s^2, gravitational acceleration of the earth
f = 2.1e6;          % N, thrust of the first stage of Titan II Rocket
hscale = 8440 ;    % m, atmospheric scale height
beta = h/hscale;    % nondim, constant used to reduce EOM equations
rhoref = 1.225 ;    % kg/m^3, reference density
A = 7.069;          % m^2, aerodynamic reference area
rad2deg = 180/pi;

%% Boundary Conditions
global xbar0 ybar0 Vxbar0 Vybar0 ybarf Vxbarf Vybarf     % pass to BC function

% initial conditions`
xbar0 = 0;         % meters
ybar0 = 0;
Vxbar0 = 0;        % meters/second
Vybar0 = 0;

% final conditions
ybarf = h/h;         % meters, final altitude
Vxbarf = Vc/Vc;       % m/s, final downrange velocity
Vybarf = 0;        % final vertical velocity

%% Parameters No drag, constant mass approximation
m0 = 6088;      % kg, average mass of a Titan II
CD = 0;         % Drag coefficient for a no drag case
mdot = 0;       % kg/s, mass flow rate for the constant mass
% eta is a constant combining constant terms for drag EOM
eta = rhoref*CD*A/2 ;

%% Initial Guess No drag, constant mass approximation
t0 = 0;
% list initial conditions in yinit, use 0 if unknown
yinit = [xbar0 ybar0 Vxbar0 Vybar0 0 -1 0];        % guess for initial state and costate
tf_guess = 700;                     % sec, initial guess for final time

% in order to leave tf as a free variable to optimize, we parameterize by
% tf. tau = time/tf
Nt = 80;                    % number of time points time is discretized to 
tau = linspace(0,1,Nt)';    % non-dimensional time vector

% use bvpinit to get the initial guess
% requires non-dimensional time vector; initial states, costates and 
% guesses; and guess for the final time
solinit = bvpinit(tau, yinit, tf_guess);

%% solution No drag, constant mass approximation

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
x_sol = z(1,:)*h/1000;
y_sol = z(2,:)*h/1000;
vx_sol = z(3,:)*Vc/1000;
vy_sol = z(4,:)*Vc/1000;
lambda2_bar_sol = z(5,:);
lambda3_bar_sol = z(6,:);
lambda4_bar_sol = z(6,:);

%% Parameters for Varying Masss and No Drag
m0 = 117020;
mdot = (117020 - 4760)/139;
delta_tf = 115;

%% Solution for for Varying Masss and No Drag
solinit_mass = solinit;
% update guess with solution to 7 state and costate variables from previous
% approximation
solinit_mass.y = z;
solinit_mass.parameters(1) = tf-delta_tf;

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
    global xbar0 ybar0 Vxbar0 Vybar0 ybarf Vxbarf Vybarf Vc f
    global eta beta g_accel m0 mdot
    mf = m0-abs(mdot)*tf;
    PSI = [ Y0(1) - xbar0;
            Y0(2) - ybar0;
            Y0(3) - Vxbar0;
            Y0(4) - Vybar0;
            Yf(2) - ybarf;
            Yf(3) - Vxbarf;
            Yf(4) - Vybarf;
            (-sqrt(Yf(6)^2+Yf(7)^2)*f/mf/Vc...
            -(Yf(6)*Yf(3))*eta*exp(-beta)*sqrt(Yf(3)^2)*Vc/mf-Yf(7)...
            *g_accel/Vc)*tf+1
          ];
end

function dX_dtau = ascent_odes_tf(tau, X, tf)
    global g_accel Vc h eta beta f m0 mdot
    m = m0-abs(mdot)*tau*tf;
    xbardot = X(3)*Vc/h;
    ybardot = X(4)*Vc/h;
    Vxbardot = (f/Vc*(-X(6)/sqrt(X(6)^2+X(7)^2))...
        -eta*exp(-X(2)*beta)*X(3)*sqrt(X(3)^2+X(4)^2)*Vc)/m;
    Vybardot = (f/Vc*(-X(7)/sqrt(X(6)^2+X(7)^2))...
        -eta*exp(-X(2)*beta)*X(4)*sqrt(X(3)^2+X(4)^2)*Vc)/m - g_accel/Vc;
    if sqrt(X(3)^2+X(4)^2) == 0
        lambda2_bar = 0;
        lambda3_bar = 0;
        lambda4_bar = -X(5)*Vc/h;
    else
        lambda2_bar = ...
            -(X(6)*X(3)+X(7)*X(4))*eta*beta*exp(-X(2)*beta)*sqrt(X(3)^2+X(4)^2)*Vc/m;
        lambda3_bar = eta*exp(-X(2)*beta)*Vc*(X(6)*(2*X(3)^2+X(4)^2)...
            + X(7)*X(3)*X(4))/sqrt(X(3)^2+X(4)^2)/m;
        lambda4_bar = -X(5)*Vc/h + eta*exp(-X(2)*beta)*Vc*(X(7)*(X(3)^2 ...
            + 2*X(4)^2)...
            + X(6)*X(3)*X(4))/sqrt(X(3)^2+X(4)^2)/m;
    end
    dX_dtau = tf*[xbardot; ybardot; Vxbardot; Vybardot; lambda2_bar; lambda3_bar; lambda4_bar];
end