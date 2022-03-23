clear;clc;close all;
syms x L F0
phi(1) = 1 - 3.*(x/L).^2 + 2.*(x/L).^2;
phi(2) = -x.*(1-x/L).^2;
phi(3) = 3.*(x/L).^2 - 2.*(x/L).^3;
phi(4) = -x.*((x/L).^2 - (x/L));

qf = int(F0*dirac(x-L/2),x,0,L);
q = sym(zeros(4,1));
qd = sym(zeros(4,1));
qt = sym(zeros(4,1));
for ii = 1:4
    q(ii) = int(F0*phi(ii),x,0,L);
    qd(ii) = int(F0*dirac(x-L/2)*phi(ii),x,0,L);
    qt(ii) = int(qf*phi(ii),x,0,L);
end
disp(q)
disp(qd)
disp(qt)

