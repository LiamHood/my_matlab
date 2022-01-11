%% Example 4.3
clear; clc; close all;
% minimize J = 1/2 * int(x^2 + P*u^2,0,1)
    % x^2 is square of error
    % u^2 is control effort
    % xdot = u
    % x(0) = 1
% Let P = 1
% Form Hamiltonian
    % H = (x^2)/2 + (u^2)/2 + lambda*u
% Write Euler-Lagrange equations
    % lambdadot = -Hx
    % Hu = 0
    % lambdadot = -delH/delx = -diff(H,x) = -x
    % delH/delu = diff(H,u) = lambda + uds = 0
    % u = -lambda
% differential form of transversability condition
    % Hf*dtf - lambdaf'*dxf + dphif = 0
    % dpsi = 0
    % THIS problem 
        % phi = 0
        % psi = tf - 1 = 0
    % dpsi = dtf = 0
    % lambdaf = lambda(tf) = delphi/delx(tf) = 0
% well defined TPBVP
    % xdot = -lambda
    % lambdadot = -x
    % t0 = 0
    % x0 = 1
    % tf = 1
    % lambdaf = 0
% state-variable form 
    % z = [z1; z2] = [x; lambda]
    % zdot = [z1dot; z2dot] = [-z2; -z1]
    % z1(0) = 1
    % z2(1) = 0
% analytical solution
    % z1dd = -z2dot = z1
    % z1dd - z1 = 0
    % z1 = c1*exp(t) + c2*exp(-t)
    % z2 = -z1dot = -c1*exp(t) + c2*exp(-t)
    % z1(0) = c1 + c2 = 1
    % z2(1) = -c1*exp(1) + c2*exp(-1) = 0
    % c1 = 1/(1+exp(2))
    % c2 = exp(2)/(1 + exp(2))

    
syms x u lambda
H = (x^2)/2 + (u^2)/2 + lambda*u;
lambdadot = -diff(H,x);
Hu = diff(H,u);
u = solve(Hu == 0,u);
