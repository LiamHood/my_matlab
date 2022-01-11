clear; close all; clc;
syms x xa xb u1 u2 u3 xn xp u4 u5 q11 u6 u7 u8 u9

k = @(psii,psij,xa_val,xb_val) int(4*diff(psii)*diff(psij)+2*psii*psij,x,xa_val,xb_val);
f = @(psii) int(psii*(x^2+1),x,xa,xb);

psi1a = (xb-x)/(xb-xa);
psi2a = (x-xa)/(xb-xa);

Ka = [  k(psi1a,psi1a,xa,xb), k(psi1a,psi2a,xa,xb);
        k(psi2a,psi1a,xa,xb), k(psi2a,psi2a,xa,xb)];
Fa = [f(psi1a);f(psi2a)];
Ua = [u1;u2];
ca = inv([1,xa;1,xb])*Ua;
Qa_def = [-4*diff(ca(1)+ca(2)*x);4*diff(ca(1)+ca(2)*x)];
Qa_equil = Ka*Ua-Fa;

disp(Ka);
disp(Fa);
disp(Qa_equil);
disp(Qa_def);

xbar = x-xa;
h = xb-xa;
psib(1) = (1-xbar/h)*(1-2*xbar/h);
psib(2) = 4*xbar/h*(1-xbar/h);
psib(3) = -xbar/h*(1-2*xbar/h);

Kb = sym(zeros(3));
nodesb = [xa, xa+(xa+xb)/2, xb];
for ii = 1:3
    for jj = 1:3
        Kb(ii,jj) = subs(k(psib(ii),psib(jj),xa,xb),x,nodesb(ii));
    end
end

Fb = [subs(f(psib(1)),x,xa);subs(f(psib(2)),x,xa+(xa+xb)/2);subs(f(psib(3)),x,xb)];
Ub = [u1;u2;u3];
Qb_equil = Kb*Ub-Fb;
cb = inv([1,xa,xa^2;1,xa+(xb+xa)/2,(xa+(xb+xa)/2)^2;1,xb,xb^2])*Ub;
Qb_def = [subs(-4*diff(cb(1)+cb(2)*x+cb(3)*x^2),x,xa);
    subs(-4*diff(cb(1)+cb(2)*x+cb(3)*x^2),x,xp)+subs(4*diff(cb(1)+cb(2)*x+cb(3)*x^2),x,xn);
    subs(4*diff(cb(1)+cb(2)*x+cb(3)*x^2),x,xb)];

disp(Kb);
disp(Fb);
disp(Qb_def);
disp(Qb_equil);

K2a(:,:,1) = subs(Ka,[xa,xb],[0,2]);
K2a(:,:,2) = subs(Ka,[xa,xb],[2,4]);
K2a(:,:,3) = subs(Ka,[xa,xb],[4,6]);
K2a(:,:,4) = subs(Ka,[xa,xb],[6,8]);
U2a = [0;u2;u3;u4;u5];
F2a(1,:) = subs(Fa,[xa,xb],[0,2]);
F2a(2,:) = subs(Fa,[xa,xb],[2,4]);
F2a(3,:) = subs(Fa,[xa,xb],[4,6]);
F2a(4,:) = subs(Fa,[xa,xb],[6,8]); 
K2ag = zeros(5);
F2ag = zeros(5,1);
nodesa = [0,2,4,6,8];
for ee = 1:4
    for ii = 1:2
        for jj = 1:2
            K2ag(ee-1+ii,ee-1+jj) = K2ag(ee-1+ii,ee-1+jj) + K2a(ii,jj,ee);
        end
        F2ag(ee-1+ii) = subs(F2a(ee,ii),x,nodesa(ee-1+ii));
    end
end
Q2a = [q11;0;0;0;200];
ans2a = solve(K2ag*U2a==F2ag+Q2a,[q11,u2,u3,u4,u5]);
U2asol = [0;ans2a.u2;ans2a.u3;ans2a.u4;ans2a.u5];

K2b(:,:,1) = subs(Kb,[xa,xb],[0,2]);
K2b(:,:,2) = subs(Kb,[xa,xb],[2,4]);
K2b(:,:,3) = subs(Kb,[xa,xb],[4,6]);
K2b(:,:,4) = subs(Kb,[xa,xb],[6,8]);
U2b = [0;u2;u3;u4;u5;u6;u7;u8;u9];
F2b(1,:) = subs(Fb,[xa,xb],[0,2]);
F2b(2,:) = subs(Fb,[xa,xb],[2,4]);
F2b(3,:) = subs(Fb,[xa,xb],[4,6]);
F2b(4,:) = subs(Fb,[xa,xb],[6,8]); 
K2bg = zeros(9);
F2bg = zeros(9,1);
nodesb = [0,1,2,3,4,5,6,7,8];
for ee = 1:4
    for ii = 1:3
        efac = 2*(ee-1);
        for jj = 1:3
            K2bg(efac+ii,efac+jj) = K2bg(efac+ii,efac+jj) + K2b(ii,jj,ee);
        end
        F2bg(efac+ii) = subs(F2b(ee,ii),x,nodesb(efac+ii));
    end
end
Q2b = [q11;0;0;0;0;0;0;0;200];
ans2b = solve(K2bg*U2b==F2bg+Q2b,[q11,u2,u3,u4,u5,u6,u7,u8,u9]);
U2bsol = [0;ans2b.u2;ans2b.u3;ans2b.u4;ans2b.u5;ans2b.u6;ans2b.u7;ans2b.u8;ans2b.u9];

disp("Q11 equil")
disp(vpa(ans2a.q11))
disp("Q42 equil")
disp(200)
disp("Q11 def")
disp(vpa(-4*U2asol(1)*subs(diff(subs(psi1a,[xa,xb],[0,2]),x),x,0)))
disp("Q42 def")
disp(vpa(4*U2asol(5)*subs(diff(subs(psi2a,[xa,xb],[6,8]),x),x,8)))

disp("Q11 equil")
disp(vpa(ans2b.q11))
disp("Q43 equil")
disp(200)
disp("Q11 def")
disp(vpa(-4*U2bsol(1)*subs(diff(subs(psib(1),[xa,xb],[0,2]),x),x,0)))
disp("Q43 def")
disp(vpa(4*U2bsol(5)*subs(diff(subs(psib(3),[xa,xb],[6,8]),x),x,8)))

clear x
u_exact = @(x) .2074696334*exp(.7071067810*x)-2.707469632*exp(-.7071067810*x)+.5*x.^2+2.5;
x = linspace(0,8,41);
xele = [linspace(0,2,11);linspace(2,4,11);linspace(4,6,11);linspace(6,8,11)];
U2acont = zeros(41,1);

for ee = 1:4
    upper = (10*(ee))+1;
    lower = (10*(ee-1))+1;
    U2acont(lower:upper) = U2asol(ee)*subs(subs(psi1a,[xa,xb],[nodesa(ee),nodesa(ee+1)]),xele(ee,:)) + ...
                             U2asol(ee+1)*subs(subs(psi2a,[xa,xb],[nodesa(ee),nodesa(ee+1)]),xele(ee,:));
    U2bcont(lower:upper) = U2bsol(ee)*subs(subs(psib(1),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:)) + ...
                           U2bsol(ee+1)*subs(subs(psib(2),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:)) + ...
                           U2bsol(ee+2)*subs(subs(psib(2),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:));
end
figure
hold on
plot(x, u_exact(x))
plot(x, U2acont)
plot(x, U2bcont)
