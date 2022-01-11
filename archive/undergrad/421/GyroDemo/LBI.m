function L_B_I = LBI(phi, theta, psi, sequence)

% 3-1-3 sequence
if sequence == 313
    L_B_I = 1/sin(theta)*[sin(psi) cos(psi) 0;
        cos(psi)*sin(theta) -sin(psi)*sin(theta) 0;
        -sin(psi)*cos(theta) -cos(psi)*cos(theta) sin(theta)];
elseif sequence == 321
    % 3-2-1 sequence
    L_B_I = 1/cos(theta)*[cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta);
        0 cos(phi)*cos(theta) -sin(phi)*cos(theta);
        0 sin(phi) cos(phi)];
end