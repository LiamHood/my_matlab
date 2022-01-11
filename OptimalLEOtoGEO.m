% LEO to GEO Circular Orbit Radius Transfer Problem with MATLAB’s bvp4c
%
% by Peter J. Edelman and Kshitij Mall
%
% This script uses MATLAB’s bvp4c to solve the problem of finding the
% optimal transfer trajectory for launch from a 300 kilometer altitude LEO
% to a 35,786 kilometer altitude GEO orbit.
% A linearly time-varying mass is assumed.
% In addition to this script, we use two functions:
% one to provide the differential equations and another that gives the
% boundary conditions:
%
% This file: Max_Xfer.m
% State and Costate Equations: transfer_odes.m
% Boundary Conditions: transfer_bcs.m
%
% Clear figure, command and workspace histories
clear all; close all; clc;
% Define parameters of the problem
% pass these parameters to the DE and BC functions
global eta mdot Tbar tf m0
h0 = 300; % Initial altitude [km]
rearth = 6378; % Mean radius of the Earth [km]
mu = 3.986004418e5; % Gravitational parameter of Earth [kmˆ3/sˆ2]
T = 20; % Thrust of spacecraft [N]
Isp = 6000; % Specific Impulse [sec]
mdot = T/Isp/9.80665; % Mass flow rate [kg/s]
m0 = 1500; % Initial spacecraft mass [kg]

%---------------------------------------------------------------------------
%% Boundary Conditions
%---------------------------------------------------------------------------
global r0bar u0bar v0bar theta0 ufbar % pass these BCs to the BC function
% Initial conditions
% Start from 300 km altitude circular orbit
% All Boundary Conditions are nondimensional
r0 = rearth+h0; % Initial circular orbit radius [km]
u0 = 0; % Initial radial component of velocity [km/s]
v0 = sqrt(mu/r0); % Initial circular orbit speed [km/s]
theta0 = 0; % Initial longitude [rad]
t0 = 0; % Initial time [sec]
% Final conditions
uf = 0; % Final radial component of velocity [km/s]
days = .1; % Initial amount of days to thrust. (0.1 days)
tf = 24*3600*days; % Days in seconds

%---------------------------------------------------------------------------
%% Nondimensionalize Boundary Conditions and Important Parameters
%---------------------------------------------------------------------------
% Scaling constant eta
eta = v0*tf/r0;
% Scaled thrust Tbar
Tbar = T*tf/(v0*1000);
% Scaled Initial Conditions
r0bar = r0/r0;
u0bar = u0/v0;
v0bar = v0/v0;
% Scaled Final Conditions
ufbar = uf/v0;

%---------------------------------------------------------------------------
%% Solve the TPBVP
%---------------------------------------------------------------------------
% Nt is the number of points that the TPBVP will be discretized into. First
% 1000 mesh points will be used at each iteration to reduce computational
% time. After the very last iteration, the solution will be spline
% interpolated to 5000 mesh points and solved again, so as to ensure smooth
% plots.
Nt = 1000;
tau = linspace(0,1,Nt); % nondimensional time vector
% list initial conditions in yinit
yinit = [r0bar u0bar v0bar 0 0 -1 0];
% Create an initial guess of the solution using the MATLAB function
% bvpinit, which requires as inputs the (nondimensional) time vector,
% initial states (or guesses, as applicable).
% The initial guess of the solution is stored in the structure solinit.
solinit = bvpinit(tau,yinit);
% Call bvp4c to solve the TPBVP. Point the solver to the functions
% containing the differential equations and the boundary conditions and
% provide it with the initial guess of the solution.
sol = bvp4c(@transfer_odes,@transfer_bcs,solinit);
% Evaluate the solution at all times in the nondimensional time vector tau
% and store the state variables in the matrix Z.
Z = deval(sol,tau);
% Implement for-loop to increment final time by 0.1 days each iteration
for days = .2:.1:2
    % update final time
    tf = 24*3600*days;
    % update Tbar and eta, which are functions of tf
    Tbar = T*tf/(v0*1000);
    eta = v0*tf/r0;
    % store the previous solution as the guess for the next iteration
    solinit.y = Z;
    solinit.x = tau;
    % solve the TPBVP
    sol = bvp4c(@transfer_odes,@transfer_bcs,solinit);
    % store the newly found solution in a matrix ’Z’
    Z = deval(sol,tau);
end
% for-loop to increment final time by 0.05 days each iteration for days 2 -
% 3.95. Due to numerical sensitivity, the step size had to be reduced.
for days = 2.05:.05:3.95
    % update final time
    tf = 24*3600*days;
    % update Tbar and eta, which are functions of tf
    Tbar = T*tf/(v0*1000);
    eta = v0*tf/r0;
    % store the previous solution as the guess for the next iteration
    solinit.y = Z;
    solinit.x = tau;
    % solve the TPBVP
    sol = bvp4c(@transfer_odes,@transfer_bcs,solinit);
    % store the newly found solution in a matrix ’Z’
    Z = deval(sol,tau);
end
% Final iteration to get to GEO (still 1000 mesh points)
days = 3.97152;
% update final time
tf = 24*3600*days;
% update Tbar and eta, which are functions of tf
Tbar = T*tf/(v0*1000);
eta = v0*tf/r0;
% store the previous solution as the guess for the next iteration
solinit.y = Z;
solinit.x = tau;
% solve the TPBVP
sol = bvp4c(@transfer_odes,@transfer_bcs,solinit);
% store the newly found solution in a matrix ’Z’
Z = deval(sol,tau);

%% Interpolate with spline to get 5000 mesh points
% Create 5000 point time vector
tau2 = linspace(0,1,5*Nt);
% Interpolate solution with spline using the new time vector
for q = 1:7
    Z2(q,:) = spline(tau,Z(q,:),tau2);
end
% Store the new solution in tau and Z
tau = tau2;
Z = Z2;

%% Solve TPBVP again with the 5000 mesh points for "smooth" answer
solinit.y = Z;
solinit.x = tau;
sol = bvp4c(@transfer_odes,@transfer_bcs,solinit);
% store the newly found solution in a matrix ’Z’
Z = deval(sol,tau);
% Convert back to dimensional time for plotting
time = t0+tau*(tf-t0);
% Extract the solution for each state variable from the matrix Z and
% convert them back into dimensional units by multiplying each by their
% respective scaling constants.
r_sol = Z(1,:)*r0; % Radius [km]
u_sol = Z(2,:)*v0; % Radial component of velocity [km/s]
v_sol = Z(3,:)*v0; % Tangential component of velocity [km/s]
theta_sol = Z(4,:); % Angle between x-axis and radius vector [rad]
lambda_rbar_sol = Z(5,:); % 1st costate
lambda_ubar_sol = Z(6,:); % 2nd costate
lambda_vbar_sol = Z(7,:); % 3rd costate
% Displays final value of orbit radius
final_radius = r_sol(end);

%---------------------------------------------------------------------------
%% Plots
%---------------------------------------------------------------------------
figure(1)
plot(r_sol.*cos(theta_sol), r_sol.*sin(theta_sol),'k')
xlabel('x-direction [km]')
ylabel('y-direction [km]')
axis equal
% Plot the steering angle as a function of time
figure(2)
plot(time/(3600*24), atand(lambda_ubar_sol./lambda_vbar_sol),'k')
xlabel('Time [days]')
ylabel('\theta(t) [deg]')
xlim([time(1) time(end)]/(3600*24))
function dX_dtau = transfer_odes(tau,x)
%
% State and Costate Differential Equation Function for the Maximal Radius
% Low-Thrust Orbit Transfer
%
%
% The independent variable here is the nondimensional time, tau, and the
% state vector is X
% Note that the state vector X has components
% x(1) = rbar, nondimensional radius from center of Earth
% x(2) = ubar, nondimensional radial component of velocity
% x(3) = vbar, nondimensional tangential component of velocity
% x(4) = theta, angle between inertial x-axis and radius vector
% x(5) = lambda_r_bar, first costate
% x(6) = lambda_u_bar, second costate
% x(7) = lambda_v_bar, third costate
% pass in values of relevant constants as global variables
global eta mdot Tbar tf m0
m = m0-abs(mdot)*tau*tf;
drbar_dtau = x(2)*eta;
dubar_dtau = x(3)^2/x(1)*eta-eta/x(1)^2 ...
-Tbar/m*(x(6)/sqrt(x(6)^2+x(7)^2));
dvbar_dtau = -x(2)*x(3)/x(1)*eta-Tbar/m*(x(7)/sqrt(x(6)^2+x(7)^2));
dtheta_dtau = x(3)/x(1)*eta;
dlambda_r_bar_dtau = x(6)*(x(3)^2/x(1)^2*eta-2*eta/x(1)^3) ...
-x(7)*x(2)*x(3)/x(1)^2*eta;
dlambda_u_bar_dtau = -x(5)*eta+x(7)*x(3)/x(1)*eta;
dlambda_v_bar_dtau = -x(6)*2*x(3)/x(1)*eta+x(7)*x(2)/x(1)*eta;
dX_dtau = [drbar_dtau; dubar_dtau; dvbar_dtau; dtheta_dtau; ...
dlambda_r_bar_dtau; dlambda_u_bar_dtau; dlambda_v_bar_dtau];
end

function PSI = transfer_bcs(Y0,Yf)
%
% Boundary Condition Function for the Maximal Radius Low-Thrust Orbit
% Transfer
%
    % pass in values of boundary conditions and other constants as global variables
    global r0bar u0bar v0bar ufbar theta0
    % Create a column vector with the 7 boundary conditions on the state &
    % costate variables
    PSI = [
        Y0(1)-r0bar % Initial Condition
        Y0(2)-u0bar % Initial Condition
        Y0(3)-v0bar % Initial Condition
        Y0(4)-theta0 % Initial Condition
        Yf(2)-ufbar % Final Condition
        Yf(3)-sqrt(1/Yf(1)) % Final Condition
        -Yf(5)+1/2*Yf(7)/Yf(1)^(3/2)-1
        ];% Final Condition
end