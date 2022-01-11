clear; close all; clc;
syms x xa xb u1 u2 u3 xn xp u4 u5 q11 u6 u7 u8 u9
disp("Problem 1")
k = @(psii,psij,xa_val,xb_val) int(4*diff(psii)*diff(psij)+2*psii*psij,x,xa_val,xb_val);
f = @(psii) int(psii*(x^2+1),x,xa,xb);

xbar = x-xa;
h = xb-xa;

psia(1) = (xb-x)/(h);
psia(2) = (xbar)/(h);

psib(1) = (1-xbar/h)*(1-2*xbar/h);
psib(2) = 4*xbar/h*(1-xbar/h);
psib(3) = -xbar/h*(1-2*xbar/h);

Ka = [  k(psia(1),psia(1),xa,xb), k(psia(1),psia(2),xa,xb);
        k(psia(2),psia(1),xa,xb), k(psia(2),psia(2),xa,xb)];
Fa = [f(psia(1));f(psia(2))];
Ua = [u1;u2];
Qa_equil = Ka*Ua-Fa;
Qa_def = [(u1-u2)/h;(u2-u1)/h];
disp("Stiffness Matrix")
disp(Ka)
disp("Source Vector")
disp(Fa)
disp("Secondary Nodal Degrees of Freedom")
disp("By Definition")
disp(Qa_def)
disp("By Equilibrium")
disp(Qa_equil)

Kb = sym(zeros(3));
nodesb = [xa, xa+(xa+xb)/2, xb];
for ii = 1:3
    for jj = 1:3
%         Kb(ii,jj) = subs(k(psib(ii),psib(jj),xa,xb),x,nodesb(ii));
        Kb(ii,jj) = k(psib(ii),psib(jj),xa,xb);
    end
end
% Fb = [subs(f(psib(1)),x,nodesb(1));subs(f(psib(2)),x,nodesb(2));subs(f(psib(3)),x,nodesb(3))];
Fb = [subs(f(psib(1)),x,nodesb(1));subs(f(psib(2)),x,nodesb(2));subs(f(psib(3)),x,nodesb(3))];
Ub = [u1;u2;u3];
Qb_equil = Kb*Ub-Fb;
Qb_def = [8*(u1-u2)/h;0;8*(u3-u2)/h];
disp("Stiffness Matrix")
disp(Kb)
disp("Source Vector")
disp(Fb)
disp("Secondary Nodal Degrees of Freedom")
disp("By Definition")
disp(Qb_def)
disp("By Equilibrium")
disp(Qb_equil)

disp("Problem 2")
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
%         F2ag(ee-1+ii) = subs(F2a(ee,ii),x,nodesa(ee-1+ii));
        F2ag(ee-1+ii) = F2a(ee,ii);
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
        F2bg(efac+ii) = F2b(ee,ii);
    end
end
Q2b = [q11;0;0;0;0;0;0;0;200];
ans2b = solve(K2bg*U2b==F2bg+Q2b,[q11,u2,u3,u4,u5,u6,u7,u8,u9]);
U2bsol = [0;ans2b.u2;ans2b.u3;ans2b.u4;ans2b.u5;ans2b.u6;ans2b.u7;ans2b.u8;ans2b.u9];
disp("Linear Elements")
disp("Q11 equil")
disp(vpa(ans2a.q11))
disp("Q42 equil")
disp(200)
disp("Q11 def")
disp(vpa(-4*U2asol(1)*subs(diff(subs(psia(1),[xa,xb],[0,2]),x),x,0)))
% disp(vpa(subs(Qa_def(1),[u1,u2,xa,xb],[U2asol(1),U2asol(2),0,2])))
disp("Q42 def")
disp(vpa(4*U2asol(5)*subs(diff(subs(psia(2),[xa,xb],[6,8]),x),x,8)))
% disp(vpa(subs(Qa_def(2),[u1,u2,xa,xb],[U2asol(4),U2asol(5),6,8])))

disp("Quadratic Elements")
disp("Q11 equil")
disp(vpa(ans2b.q11))
disp("Q43 equil")
disp(200)
disp("Q11 def")
disp(vpa(-4*U2bsol(1)*subs(diff(subs(psib(1),[xa,xb],[0,2]),x),x,0)))
% disp(vpa(subs(Qb_def(1),[u1,u2,xa,xb],[U2bsol(1),U2bsol(2),0,2])))
disp("Q43 def")
disp(vpa(4*U2bsol(5)*subs(diff(subs(psib(3),[xa,xb],[6,8]),x),x,8)))
% disp(vpa(subs(Qb_def(3),[u2,u3,xa,xb],[U2bsol(7),U2bsol(8),6,8])))

clear x
u_exact = @(x) .2074696334*exp(.7071067810*x)-2.707469632*exp(-.7071067810*x)+.5*x.^2+2.5;
x = linspace(0,8,41);
xele = [linspace(0,2,11);linspace(2,4,11);linspace(4,6,11);linspace(6,8,11)];
U2acont = zeros(41,1);

for ee = 1:4
    upper = (10*(ee))+1;
    lower = (10*(ee-1))+1;
    U2acont(lower:upper) = U2asol(ee)*subs(subs(psia(1),[xa,xb],[nodesa(ee),nodesa(ee+1)]),xele(ee,:)) + ...
                             U2asol(ee+1)*subs(subs(psia(2),[xa,xb],[nodesa(ee),nodesa(ee+1)]),xele(ee,:));
    U2bcont(lower:upper) = U2bsol(2*ee-1)*subs(subs(psib(1),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:)) + ...
                           U2bsol(2*ee)*subs(subs(psib(2),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:)) + ...
                           U2bsol(2*ee+1)*subs(subs(psib(3),[xa,xb],[nodesb(2*ee-1),nodesb(2*ee+1)]),xele(ee,:));
end
figure
hold on
plot(x, u_exact(x))
plot(x, U2acont)
plot(x, U2bcont)
legend("Exact","Linear","Quadratic")

%% Problem 3
clear;

disp("Problem 3")

syms xbar xa xb u1 u2 u3 u4
Es = 30e6;
Ea = 10e6;
d = [2.5,2,1.5];
A = (d./2).^2.*pi;
h = [10,6,10];
k = @(h,E,A,psii,psij) int(E*A*diff(psii(h))*diff(psij(h)),xbar,0,h);
psi{1} = @(h) 1-xbar/h;
psi{2} = @(h) xbar/h;
K{1}= zeros(2);
for ii = 1:2
    for jj = 1:2
        K{1}(ii,jj) = k(h(1),Es,A(1),psi{ii},psi{jj});
    end
end
K{2} = zeros(2);
for ii = 1:2
    for jj = 1:2
        K{2}(ii,jj) = k(h(2),Ea,A(2),psi{ii},psi{jj});
    end
end
K{3} = zeros(2);
for ii = 1:2
    for jj = 1:2
        K{3}(ii,jj) = k(h(3),Es,A(3),psi{ii},psi{jj});
    end
end
disp("K element 1")
disp(K{1})
disp("K element 2")
disp(K{2})
disp("K element 3")
disp(K{3})
Kg = zeros(4);
for ee = 1:3
    for ii = 1:2
        for jj = 1:2
            ofs = ee - 1;
            Kg(ofs+ii,ofs+jj) = Kg(ofs+ii,ofs+jj) + K{ee}(ii,jj);
        end
    end
end
U = [u1;u2;u3;u4];
Q = [-100;0;-200;300];
disp("Global K")
disp(Kg)
disp("Global U")
disp(U)
disp("Global Q")
disp(Q)
disp("Equation")
disp("KU=Q")
disp(Kg*U==Q)
disp("U")
disp(pinv(Kg)*Q)

%% Problem 4
clear;

disp("Problem 4")

syms xbar xa xb u1 u2 u3 u4 Q1 Q2
h = pi/3;
k = @(h,psii,psij) int(diff(psii(h))*diff(psij(h)),xbar,0,h);
f = @(h,psii) int(psii(h)*(sin(pi*xbar/2)),xbar,0,h);
psi{1} = @(h) 1-xbar/h;
psi{2} = @(h) xbar/h;
K{1}= zeros(2);
F{1}=zeros(2,1);
for ii = 1:2
    for jj = 1:2
        K{1}(ii,jj) = k(h,psi{ii},psi{jj});
    end
    F{1}(ii) = f(h,psi{ii});
end
K{2} = zeros(2);
F{2}=zeros(2,1);
for ii = 1:2
    for jj = 1:2
        K{2}(ii,jj) = k(h,psi{ii},psi{jj});
    end
    F{2}(ii) = f(h,psi{ii});
end
K{3} = zeros(2);
F{3}=zeros(2,1);
for ii = 1:2
    for jj = 1:2
        K{3}(ii,jj) = k(h,psi{ii},psi{jj});
    end
    F{2}(ii) = f(h,psi{ii});
end
disp("K element 1")
disp(K{1})
disp("K element 2")
disp(K{2})
disp("K element 3")
disp(K{3})
Kg = zeros(4);
Fg = zeros(4,1);
for ee = 1:3
    ofs = ee - 1;
    for ii = 1:2
        for jj = 1:2
            Kg(ofs+ii,ofs+jj) = Kg(ofs+ii,ofs+jj) + K{ee}(ii,jj);
        end
        Fg(ofs+ii) = Fg(ofs+ii) + F{ee}(ii);
    end
end
U = [0;u2;u3;0];
Q = [Q1;0;0;Q2];
q_equil = [Kg(1,1),Kg(1,4);Kg(4,1),Kg(4,4)]*[0;0]-[Fg(1);Fg(4)];
u_equil = inv(Kg(2:3,2:3))*(Fg(2:3)+Q(2:3));
U = [0;u_equil(1);u_equil(2);0];
Q = [q_equil(1);0;0;q_equil(2)];
disp("Global K")
disp(Kg)
disp("Global U")
disp(U)
disp("Global Q")
disp(Q)

clear x
u_exact = @(x) (4*sin(.5*pi*x)/pi^2)-(4*sin(.5*pi^2)*x/pi^3);
x = linspace(0,pi,31);
xele = [linspace(0,h,11);linspace(0,h,11);linspace(0,h,11)];
Ucont = zeros(31,1);
for ee = 1:3
    upper = (10*(ee))+1;
    lower = (10*(ee-1))+1;
    Ucont(lower:upper) = U(ee)*subs(psi{1}(h),xbar,xele(ee,:)) +...
                           U(ee+1)*subs(psi{2}(h),xbar,xele(ee,:));
end
figure
hold on
plot(x, u_exact(x))
plot(x, Ucont)
legend("Exact","Linear")

%% Problem 5
clear;
disp("Problem 5")

syms xbar xa xb u1 u2 u3 u4 u5 Q1 Q2 h1 h2 h3 h4 h5 h6
h = [h1;h2;h3;h4;h5;h6];
ke = [10,20,30,40,50,60];
k = @(h,Ke,psii,psij) int(Ke*diff(psii(h))*diff(psij(h)),xbar,0,h);
psi{1} = @(h) 1-xbar/h;
psi{2} = @(h) xbar/h;
K{1}= sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{1}(ii,jj) = k(h(1),ke(1),psi{ii},psi{jj});
    end
end
K{2} = sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{2}(ii,jj) = k(h(2),ke(2),psi{ii},psi{jj});
    end
end
K{3} = sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{3}(ii,jj) = k(h(3),ke(3),psi{ii},psi{jj});
    end
end
K{4} = sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{4}(ii,jj) = k(h(4),ke(4),psi{ii},psi{jj});
    end
end
K{5} = sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{5}(ii,jj) = k(h(5),ke(5),psi{ii},psi{jj});
    end
end
K{6} = sym(zeros(2));
for ii = 1:2
    for jj = 1:2
        K{6}(ii,jj) = k(h(6),ke(6),psi{ii},psi{jj});
    end
end
Kg = sym(zeros(5));
Kg(1,1) = K{1}(1,1);
Kg(1,2) = K{1}(1,2);
Kg(2,1) = K{1}(2,1);
Kg(2,2) = K{1}(2,2) + K{3}(1,1) + K{2}(1,1) + K{4}(1,1);
Kg(2,3) = K{3}(1,2);
Kg(3,2) = K{3}(2,1);
Kg(2,4) = K{2}(1,2);
Kg(4,2) = K{2}(2,1);
Kg(2,5) = K{4}(1,2);
Kg(5,2) = K{4}(2,1);
Kg(3,3) = K{3}(2,2) + K{5}(1,1);
Kg(3,4) = K{5}(1,2);
Kg(4,3) = K{5}(2,1);
Kg(4,4) = K{5}(2,2) + K{6}(1,1) + K{2}(2,2);
Kg(4,5) = K{6}(1,2);
Kg(5,4) = K{6}(2,1);
Kg(5,5) = K{6}(2,2) + K{4}(2,2);

U = [0;u2;u3;u4;u5];
Q = [Q1;0;0;0;50];
disp("Global K")
disp(Kg)
disp("Global U")
disp(U)
disp("Global Q")
disp(Q)
disp("Equations")
disp("KU=Q")
disp(Kg*U==Q)
sol5 = solve(Kg*U==Q,[Q1,u2,u3,u4,u5]);
disp("Q1")
disp(sol5.Q1)
disp("U2")
disp(sol5.u2)
disp("U3")
disp(sol5.u3)
disp("U4")
disp(sol5.u4)
disp("U5")
disp(sol5.u5)

%% Problem 6
clear;
disp("Problem 6")

syms xbar w2 w3 Q1 Q2 L EI
h = L/2;
k = @(h,psii,psij) int(EI*diff(psii(h))*diff(psij(h)),xbar,0,h);
f = @(h,psii) -int(psii(h)*((xbar/2)*sin(pi*xbar/L)),xbar,0,h);
psi{1} = @(h) 1-xbar/h;
psi{2} = @(h) xbar/h;
K{1}= sym(zeros(2));
F{1}=sym(zeros(2,1));
for ii = 1:2
    for jj = 1:2
        K{1}(ii,jj) = k(h,psi{ii},psi{jj});
    end
    F{1}(ii) = f(h,psi{ii});
end
K{2} = sym(zeros(2));
F{2}=sym(zeros(2,1));
for ii = 1:2
    for jj = 1:2
        K{2}(ii,jj) = k(h,psi{ii},psi{jj});
    end
    F{2}(ii) = f(h,psi{ii});
end
Kg = sym(zeros(3));
Fg = sym(zeros(3,1));
for ee = 1:2
    ofs = ee - 1;
    for ii = 1:2
        for jj = 1:2
            Kg(ofs+ii,ofs+jj) = Kg(ofs+ii,ofs+jj) + K{ee}(ii,jj);
        end
        Fg(ofs+ii) = Fg(ofs+ii) + F{ee}(ii);
    end
end
W = [0;w2;w3] ;
Q = [Q1;0;0];
W_equil = pinv(Kg)*(Fg+W);
sol6 = solve(Kg*W==Fg+Q,[Q1,w2,w3]);
disp("Max Deflection")
disp(sol6.w3)

%% Problem 7
clear;
syms xbar u1 v1 u2 v2 u3 v3 F1x F1y F3y
disp("Problem 7")

E = 10e6;
A = [40,50,60];
h = [sqrt(20^2+12^2),sqrt(20^2+20^2),32];
theta = [atan(20/12),-atan(20/20),0];

K = @(E,A,h,theta) E*A/h * [cos(theta)^2, .5*sin(2*theta), -cos(theta)^2, -.5*sin(2*theta);
                            .5*sin(2*theta), sin(theta)^2, -.5*sin(2*theta), -sin(theta)^2;
                            -cos(theta)^2, -.5*sin(2*theta), cos(theta)^2, .5*sin(2*theta);
                            -.5*sin(2*theta), -sin(theta)^2, .5*sin(2*theta), sin(theta)^2];
Ke{1} = K(E,A(1),h(1),theta(1));
Ke{2} = K(E,A(2),h(2),theta(2));
Ke{3} = K(E,A(3),h(3),theta(3));
Kg = zeros(6);
for ii = 1:4
    for jj = 1:4
        Kg(ii,jj) = Kg(ii,jj) + Ke{1}(ii,jj);
    end
end
b2 = [1,2,5,6];
for ii = 1:4
    for jj = 1:4
        Kg(b2(ii),b2(jj)) = Kg(b2(ii),b2(jj)) + Ke{2}(ii,jj);
    end
end
b3 = [3,4,5,6];
for ii = 1:4
    for jj = 1:4
        Kg(b3(ii),b3(jj)) = Kg(b3(ii),b3(jj)) + Ke{3}(ii,jj);
    end
end
F = [F1x;F1y;-40;-30;0;F3y];
delta = [0;0;u2;v2;u3;0];
sol7 = solve(Kg*delta==F,[F1x,F1y,F3y,u2,v2,u3]);
delta = [0;0;vpa(sol7.u2);vpa(sol7.v2);vpa(sol7.u3);0];
disp("Displacement of Node 2 x")
disp(vpa(sol7.u2))
disp("Displacement of Node 2 y")
disp(vpa(sol7.v2))
disp("Reaction at Node 1 x")
disp(vpa(sol7.F1x))
disp("Reaction at Node 1  y")
disp(vpa(sol7.F1y))
disp("Reaction at Node 3  y")
disp(vpa(sol7.F3y))
delta_rel = @(theta, delta) [   cos(theta),sin(theta),0,0;
                                -sin(theta),cos(theta),0,0;
                                0,0,cos(theta),sin(theta);
                                0,0,-sin(theta),cos(theta)]*delta;
delta_rel1 = delta_rel(theta(1),delta(1:4));
delta_rel2 = delta_rel(theta(2),[delta(1:2);delta(5:6)]);
delta_rel3 = delta_rel(theta(3),delta(3:6));
disp("Stress in element 1")
disp((delta_rel1(1)+delta_rel1(3))/h(1)*E)
disp("Stress in element 2")
disp((delta_rel2(1)+delta_rel2(3))/h(2)*E)
disp("Stress in element 3")
disp((delta_rel3(1)+delta_rel3(3))/h(3)*E)
