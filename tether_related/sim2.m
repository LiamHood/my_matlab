clear; close all; clc;

%% Inputs
a0 = 6878;
e0 = .02;
i0_deg = 30;
i0 = deg2rad(i0_deg);
RAAN0_deg = 30;
RAAN0 = deg2rad(RAAN0_deg);
aop0_deg = 50;
aop0 = deg2rad(aop0_deg);
ta0_deg = 0;
ta0 = deg2rad(ta0_deg);
tend = 30;

L = 15;%km
m1 = 250;%kg
m2 = 15;%kg
mt = 20;%kg
It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
Ia = 0;
current_type = 0;
current_val = 1;
menu1 = input("Run with defaults? [Y]/N", 's');
if isempty(menu1)
    menu1 = 'Y';
end
if menu1 == 'N'
    a0 = input("Semi-Major Axis [km] (default: 6878): ");
    if isempty(a0)
        a0 = 6878;
    end
    e0 = input("Eccentricity (default: .02): ");
    if isempty(e0)
        e0 = .02;
    end
    i0_deg = input("Inclination [degree] (default: 30): ");
    if isempty(i0_deg)
        i0_deg = 30;
    end
    i0 = deg2rad(i0_deg);
    RAAN0_deg = input("Right Ascension of Ascending Node [degree] (default: 30): ");
    if isempty(RAAN0_deg)
        RAAN0_deg = 30;
    end
    RAAN0 = deg2rad(RAAN0_deg);
    aop0_deg = input("Argument of Perigee [degree] (default: 50): ");
    if isempty(aop0_deg)
        aop0_deg = 30;
    end
    aop0 = deg2rad(aop0_deg);
    ta0_deg = input("True Anomaly [degree] (default: 0): ");
    if isempty(ta0_deg)
        ta0_deg = 30;
    end
    ta0 = deg2rad(ta0_deg);
    tend = input("Time [Hours] for simulation (default: 30): ");
    if isempty(tend)
        tend = 30;
    end
    L = input("Tether length [km] (default 20): ");
    if isempty(L)
        L = 20;
    end
    m1 = input("Mass of main satellite [kg] (default: 250): ");
    if isempty(m1)
        m1 = 250;
    end
    m2 = input("Mass of subsatellite [kg] (default: 150): ");
    if isempty(m2)
        m2 = 150;
    end
    mt = input("Mass of tether [kg] (default: 20): ");
    if isempty(mt)
        mt = 20;
    end
    It = input("Moment of inertia rotating the tether (default: 8.6667e9): ");
    if isempty(It)
        It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
    end
    Ia = input("Moment of inertia about tether axis (default: 0): ");
    if isempty(Ia)
        Ia = 0;
    end
    current_type = input("Type of Current control [0: constant; 1: OML; 2: controlled] (default: 0); ");
    if isempty(current_type)
        current_type = 0;
    end
    current_val = input("Value of Current [A for constant; proportion for OML; control law] (default: 1); ");
    if isempty(current_val)
        current_val = 1;
    end
end

tspan = [0, tend*60*60];
sc_state0 = [a0; e0; i0; RAAN0; aop0; ta0];
tether_state0 = [0.1; 0.1; 0; 0; 0; 0];
tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val];

mu = 398600;
tol = 1e-12;
[ t , states] = BasicTether( tspan , sc_state0, tether_state0, ...
    tether_param, mu , tol );

%% Plot the Results
t = t/(60*60);
a = states(:,1);
e = states(:,2);
i = rad2deg(states(:,3));
RAAN = rad2deg(states(:,4));
aop = rad2deg(states(:,5));
ta = rad2deg(states(:,6));

phi = rad2deg(states(:,7));
theta = rad2deg(states(:,8));
psi = rad2deg(states(:,9));
dphi = rad2deg(states(:,10));
dtheta = rad2deg(states(:,11));

for ii = 1:length(t)
    [ rvec(:,ii) , vvec(:,ii) ] = coes2state( [sqrt(mu*a(ii)*(1-e(ii)^2)), ...
        i(ii), e(ii), RAAN(ii), aop(ii), ta(ii)] , mu );
    r(ii) = norm(rvec(:,ii));
    v(ii) = norm(vvec(:,ii));
end

figure
plot3(rvec(1,:), rvec(2,:), rvec(3,:))
xlabel("x")
ylabel("y")
zlabel("z")
axis equal

figure
subplot(2,1,1)
plot(t,r)
xlabel('time [hours]')
ylabel('Radius')

subplot(2,1,2)
plot(t,v)
xlabel('time [hours]')
ylabel('Velocity')

figure
subplot(2,1,1)
plot(t, theta)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')

subplot(2,1,2)
plot(t, phi)
xlabel('time [hours]')
ylabel('Out of Plane Libration [degree]')

figure
subplot(3,2,1)
plot(t, a)
xlabel('time [hours]')
ylabel('Semi-Major Axis [km]')

subplot(3,2,2)
plot(t, e)
xlabel('time [hours]')
ylabel('Eccentricity')

subplot(3,2,3)
plot(t, i)
xlabel('time [hours]')
ylabel('Inclination [degree]')

subplot(3,2,4)
plot(t, RAAN)
xlabel('time [hours]')
ylabel('Right Ascension of Ascending Node [degree]')

subplot(3,2,5)
plot(t, aop)
xlabel('time [hours]')
ylabel('Argument of Perigee [degree]')

subplot(3,2,6)
plot(t, ta)
xlabel('time [hours]')
ylabel('True Anomaly [degree]')



