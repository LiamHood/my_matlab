clear; clc; close all;
S = 230;
cbar = 7;
b = 34;
U1 = 170;
qbar1 = 34.3;
g = 32.1740;
n = 1;
m = 13000/g;
Ixz = 1300;
Iyy = 18800;
Izz = 47000;
ih1 = 0;

CL.alpha = 5.04;
CL.de = .40;
CL.o = 1.2;
CL.ih = .85;

Cm.alpha = -.66;
Cm.de = -.98;
Cm.o = .047;
Cm.ih = -2.1;
Cm.q = -13.5;

CD.o = .0431;
CD.alpha = 1.06;
CD.ih = 0;
CD.de = 0;

Cy.beta = -.730;
Cy.da = 0;
Cy.dr = .140;
Cy.r = .4;

Cl.beta = -.173;
Cl.da = .149;
Cl.dr = .014;
Cl.r = .450;

Cn.beta = .150;
Cn.da = -.05;
Cn.dr = -.074;
Cn.r = -.260;

r2d = 180/pi;
phiT = 0;
fprintf("1g level flight \n")
[alpha1, de1] = steadyturning2x2(CL, Cm, n, ih1, m, qbar1, S, cbar, g, U1);
fprintf("alpha: \t\t %f radian \t %f degrees \n" , alpha1, alpha1*r2d);
fprintf("elevator: \t %f radian \t %f degrees \n" , de1, de1*r2d);
thrust1 = (CD.o + CD.alpha*alpha1 + CD.ih*ih1 + CD.de*de1)*qbar1*S/cos(phiT*alpha1);
fprintf("thrust: \t %f \n", thrust1)
[beta1, da1, dr1] = lateraldirectional3x3(Cy, Cl, Cn, b, g, n, U1, Ixz, Izz, Iyy, qbar1, S);
fprintf("beta1: \t\t %f radian \t %f degrees \n" , beta1, beta1*r2d)
fprintf("aileron: \t %f radian \t %f degrees \n" , da1, da1*r2d)
fprintf("rudder: \t %f radian \t %f degrees \n" , dr1, dr1*r2d)

n = 1.5;
fprintf("\n1.5g turn \n")
[alpha2, de2] = steadyturning2x2(CL, Cm, n, ih1, m, qbar1, S, cbar, g, U1);
fprintf("alpha: \t\t %f radian \t %f degrees \n" , alpha2, alpha2*r2d);
fprintf("elevator: \t %f radian \t %f degrees \n" , de2, de2*r2d);
thrust2 = (CD.o + CD.alpha*alpha2 + CD.ih*ih1 + CD.de*de2)*qbar1*S/cos(phiT*alpha1);
fprintf("thrust: \t %f \n", thrust2)
[beta2, da2, dr2] = lateraldirectional3x3(Cy, Cl, Cn, b, g, n, U1, Ixz, Izz, Iyy, qbar1, S);
fprintf("beta1: \t\t %f radian \t %f degrees \n" , beta2, beta2*r2d)
fprintf("aileron: \t %f radian \t %f degrees \n" , da2, da2*r2d)
fprintf("rudder: \t %f radian \t %f degrees \n" , dr2, dr2*r2d)

n = 2;
fprintf("\n2g turn \n")
[alpha2, de2] = steadyturning2x2(CL, Cm, n, ih1, m, qbar1, S, cbar, g, U1);
fprintf("alpha: \t\t %f radian \t %f degrees \n" , alpha2, alpha2*r2d);
fprintf("elevator: \t %f radian \t %f degrees \n" , de2, de2*r2d);
thrust2 = (CD.o + CD.alpha*alpha2 + CD.ih*ih1 + CD.de*de2)*qbar1*S/cos(phiT*alpha1);
fprintf("thrust: \t %f \n", thrust2)
[beta2, da2, dr2] = lateraldirectional3x3(Cy, Cl, Cn, b, g, n, U1, Ixz, Izz, Iyy, qbar1, S);
fprintf("beta1: \t\t %f radian \t %f degrees \n" , beta2, beta2*r2d)
fprintf("aileron: \t %f radian \t %f degrees \n" , da2, da2*r2d)
fprintf("rudder: \t %f radian \t %f degrees \n" , dr2, dr2*r2d)


function [alpha, de] = steadyturning2x2(CL, Cm, n, ih1, m, qbar, S, cbar, g, U1)
    Cmatrix = [ CL.alpha, CL.de;
                Cm.alpha, Cm.de];
    rhs = [ n*(m*g/(qbar*S))-CL.o-CL.ih*ih1;
            -Cm.o-Cm.ih*ih1-Cm.q*(cbar*g/(2*U1^2))*(n-1/n)];
    controls = inv(Cmatrix)*rhs;
    alpha = controls(1);
    de = controls(2);
end

function [beta, da, dr] = lateraldirectional3x3(Cy, Cl, Cn, b, g, n, U1, Ixz, Izz, Iyy, qbar1, S)
    phi1 = acos(1/n);
    Cmatrix = [ Cy.beta, Cy.da, Cy.dr;
                Cl.beta, Cl.da, Cl.dr;
                Cn.beta, Cn.da, Cn.dr];
    rhs = [ -Cy.r*(b*g*sin(phi1))/(2*U1^2);
            (Izz-Iyy)*g^2*sin(phi1)^3/(qbar1*S*b*U1^2*cos(phi1))-Cl.r*b*g*sin(phi1)/(2*U1^2);
            Ixz*g^2*sin(phi1)^3/(qbar1*S*b*U1^2*cos(phi1))-Cn.r*b*g*sin(phi1)/(2*U1^2)];
    controls = inv(Cmatrix)*rhs;
    beta = controls(1);
    da = controls(2);
    dr = controls(3);
end