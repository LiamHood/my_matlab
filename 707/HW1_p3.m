fprintf("\n******Problem 3******\n")
global r2d
r2d = 180/pi;
    % 2.3-1 Find instantaneous angle of attack and angle of side slip
    % in the gusts
    
    % Given starting conditions
    v = 500;
    alpha = 8/r2d;
    beta = -5/r2d;
    % Rotation matrices
    % stability from body fixed
    C_s_bf = [cos(alpha), 0, sin(alpha);
                0, 1, 0;
                -sin(alpha), 0, cos(alpha)];
    % wind from stability
    C_w_s = [cos(beta), sin(beta), 0;
            -sin(beta), cos(beta), 0;
            0, 0, 1];
    % wind from body fixed
    C_w_bf = C_w_s*C_s_bf;
    % starting V
    vw = [v;0;0];
    vbf = C_w_bf'*vw;

    % i
    fprintf("i. Gust of 20 ft/s left to right (+velocity in y body fixed)\n")
    vbf1 = vbf + [0;20;0];
    [alpha1, beta1] = find_alpha_beta(vbf1);
    fprintf("\tAngle of attack is now %f degrees\n", alpha1)
    fprintf("\tAngle of sideslip is now %f degrees\n", beta1)
    % ii
    fprintf("ii. Gust of 50 ft/s from dead astern(+velocity in x body fixed)\n")
    vbf2 = vbf + [50;0;0];
    [alpha2, beta2] = find_alpha_beta(vbf2);
    fprintf("\tAngle of attack is now %f degrees\n", alpha2)
    fprintf("\tAngle of sideslip is now %f degrees\n", beta2)
    % iii
    fprintf("iii. Gust of 30 ft/s from right and below (-velocity in y, -velocity in z\n")
    vbf3 = vbf + [0;-30*cos(70/r2d);-30*sin(70/r2d)];
    [alpha3, beta3] = find_alpha_beta(vbf3);
    fprintf("\tAngle of attack is now %f degrees\n", alpha3)
    fprintf("\tAngle of sideslip is now %f degrees\n", beta3)

function [alpha, beta] = find_alpha_beta(vbf)
    global r2d
    syms alphas betas V
    C_s_bf = [cos(alphas), 0, sin(alphas);
                0, 1, 0;
                -sin(alphas), 0, cos(alphas)];
    % wind from stability
    C_w_s = [cos(betas), sin(betas), 0;
            -sin(betas), cos(betas), 0;
            0, 0, 1];
    % wind from body fixed
    C_w_bf = C_w_s*C_s_bf;
    sol = vpasolve([V;0;0]==C_w_bf*vbf,[alphas,betas,V]);
    alpha = sol.alphas*r2d;
    beta = sol.betas*r2d;
end