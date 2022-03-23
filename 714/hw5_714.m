clear; close all; clc;

w = 564000;
gamma = (-10:2:10)*(pi/180);
ih = 0 ;
qbar = 58;
S = 5500;
phiT = 0;
c0.d = .0751;
c0.l = .92;
c0.m = .0;
calpha.d = 1.13;
calpha.l = 5.67;
calpha.m = -1.45;
cih.d = 0;
cih.l = .75;
cih.m = -3;
cdel.d = 0;
cdel.l = .36;
cdel.m = -1.4;
dT = 5.4;
cbar = 27.3;

for ii = 1:length(gamma)
    [alpha(ii), del_e(ii), cT(ii)] = longStraight3eqn(w, gamma(ii), c0, calpha, cih, ih, cdel, qbar, S, cbar, phiT, dT);
    [alpham(ii), del_em(ii), cTm(ii)] = longStraight3by3(w, gamma(ii), c0, calpha, cih, ih, cdel, qbar, S, cbar, dT);
end
r2d = 180/pi;
disp("Problem 1")
figure
subplot(3,1,1)
hold on
plot(gamma*r2d, del_e*r2d)
plot(gamma*r2d, del_em*r2d)
hold off
xlabel('Flight Path Angle (deg)')
ylabel('Elevator Angle (deg)')
legend('With trig', 'Matrix')
subplot(3,1,2)
hold on
plot(gamma*r2d, alpha*r2d)
plot(gamma*r2d, alpham*r2d)
hold off
xlabel('Flight Path Angle (deg)')
ylabel('Angle of Attack (deg)')
legend('With trig', 'Matrix')
subplot(3,1,3)
hold on
plot(gamma*r2d, cT)
plot(gamma*r2d, cTm)
plot([-10,10],[(500/(qbar*S)),(500/(qbar*S))])
plot([-10,10],[(4*60000/(qbar*S)),(4*60000/(qbar*S))])
hold off
xlabel('Flight Path Angle (deg)')
ylabel('Thrust')
legend('With trig', 'Matrix')

thrust_valid = cTm < 4*(60000/(qbar*S)) & cTm > (500/(qbar*S)) ;
results_1_titles = ["Flight Path Angle (deg)", "Angle of Attack (deg)", "Elevator Deflection (deg)", "Coefficient of Thrust", "Is Valid"];
results_1_data = [gamma'*r2d, alpham'*r2d, del_em'*r2d, cTm', thrust_valid'];
results_1 = [results_1_titles; results_1_data];
disp(results_1)

%% Part 2
clear
U1 = 100:50:400;
gamma = 0;
phi = -5;
w = 8750;
FyT = 0;
S = 253;
qbar = .5*.002378*U1.^2;
LT = 1682;
NT = 10940;
b = 38;
DeltaND = 1000;
cy.beta = -.0105;
cy.da = 0;
cy.dr = .0021;
cl.beta = -.0029;
cl.da = .0024;
cl.dr = .0002;
cn.beta = .0018;
cn.da = .0005;
cn.dr = -.0010;

for ii = 1:length(U1)
    [beta(ii), dela(ii), delr(ii)] = OEIlatStraight3by3(w, cy, cl, cn, phi, gamma, FyT, qbar(ii), S, LT, b, NT, DeltaND);
end

disp("Problem 2")
figure
subplot(3,1,1)
hold on
plot(U1, beta)
plot([min(U1),max(U1)],[-10,-10])
plot([min(U1),max(U1)],[10,10])
xlabel('Airspeed (fps)')
ylabel('Side Slip (deg)')
hold off

subplot(3,1,2)
hold on
plot(U1, dela)
plot([min(U1),max(U1)],[-25,-25])
plot([min(U1),max(U1)],[25,25])
xlabel('Airspeed (fps)')
ylabel('Aileron Angle (deg)')
hold off

subplot(3,1,3)
hold on
plot(U1, delr)
plot([min(U1),max(U1)],[-25,-25])
plot([min(U1),max(U1)],[25,25])
xlabel('Airspeed (fps)')
ylabel('Rudder Angle (deg)')
hold off

valid = delr < 25 & delr > -25 & dela < 25 & dela > -25 & beta < 10 & beta > -10;
Vmc3by3 = min(U1(valid));
rho = .002378;
Vmc_eqn = sqrt((-2*(NT+DeltaND))/(rho*cn.dr*25*S*b));

results_2_titles = ["U (fps)", "Side Slip (deg)", "Aileron (deg)", "Rudder (deg)", "Is Valid"];
results_2_data = [U1', beta', dela', delr', valid'];
results_2 = [results_2_titles; results_2_data];
disp(results_2)

disp("Vmc from 3-by-3 solution")
disp(Vmc3by3)
disp("Vmc from 4.81")
disp(Vmc_eqn)

%% Functions
function [val_alpha, val_del_e, val_cT] = longStraight3eqn(w, gamma, c0, calpha, cih, ih, cdel, qbar, S, cbar, phiT, dT)
    syms del_s alpha_s cT_s
    eqnx = w*sin(gamma)/(qbar*S) == -(c0.d + calpha.d*alpha_s + cih.d*ih + cdel.d*del_s) + cT_s*cos(phiT + alpha_s);
    eqny = w*cos(gamma)/(qbar*S) == (c0.l + calpha.l*alpha_s + cih.l*ih + cdel.l*del_s) + cT_s*sin(phiT + alpha_s);
    eqnm = 0 == (c0.m + calpha.m*alpha_s + cih.m*ih + cdel.m*del_s)*qbar*S*cbar - cT_s*dT;
    [struct_alpha, struct_del_e, struct_cT] = vpasolve([eqnx, eqny, eqnm], [alpha_s, del_s, cT_s]);
    val_del_e = double(struct_del_e);
    val_alpha = double(struct_alpha);
    val_cT = double(struct_cT);
end

function [ alpha, del_e, cT] = longStraight3by3(w, gamma, c0, calpha, cih, ih, cdel, qbar, S, cbar, dT)
    Cmatrix = [ calpha.d, cdel.d, -1; 
                calpha.l, cdel.l, 0; 
                calpha.m, cdel.m, -dT/cbar];
    right = [   -(c0.d + cih.d*ih) - w*sin(gamma)/(qbar*S);
                -(c0.l + cih.l*ih) + w*cos(gamma)/(qbar*S);
                -(c0.m + cih.m*ih) ];
    
    
     out = (inv(Cmatrix)*right);
     alpha = out(1);
     del_e = out(2);
     cT = out(3);
end

function [beta, dela, delr] = OEIlatStraight3by3(w, cy, cl, cn, phi, gamma, FyT, qbar, S, LT, b, NT, DeltaND)
    Cmatrix = [ cy.beta, cy.da, cy.dr;
                cl.beta, cl.da, cl.dr;
                cn.beta, cn.da, cn.dr];
    right = [   -(w*sind(phi)*cosd(gamma) + FyT)/(qbar*S);
                -LT/(qbar*S*b);
                (-NT - DeltaND)/(qbar*S*b)];
    out = inv(Cmatrix)*right;
    beta = out(1);
    dela = out(2);
    delr = out(3);
end