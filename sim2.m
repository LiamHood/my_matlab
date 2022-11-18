clear; close all; clc;

tspan = [0, 1e5];

a0 = 6878;
e0 = .02;
i0 = deg2rad(30);
RAAN0 = deg2rad(30);
aop0 = deg2rad(50);
ta0 = 0;
sc_state0 = [a0; e0; i0; RAAN0; aop0; ta0];

tether_state0 = [0.1; 0.1; 0; 0; 0; 0];

L = 15;%km
m1 = 250;%kg
m2 = 15;%kg
mt = 20;%kg
It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
Ia = 0;
tether_param = [L; m1; m2; mt; It; It; Ia; 0; 1];

mu = 398600;
tol = 1e-12;
[ t , states] = BasicTether( tspan , sc_state0, tether_state0, ...
    tether_param, mu , tol );

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

figure
plot(t,r)
xlabel('time [hours]')
ylabel('Radius')

figure
plot(t,v)
xlabel('time [hours]')
ylabel('Velocity')

figure
plot(t, theta)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')

figure
plot(t, phi)
xlabel('time [hours]')
ylabel('Out of Plane Libration [degree]')

figure
plot(t, a)
xlabel('time [hours]')
ylabel('Semi-Major Axis [km]')

figure
plot(t, e)
xlabel('time [hours]')
ylabel('Eccentricity')

figure
plot(t, i)
xlabel('time [hours]')
ylabel('Inclination [degree]')

figure
plot(t, RAAN)
xlabel('time [hours]')
ylabel('Right Ascension of Ascending Node [degree]')

figure
plot(t, aop)
xlabel('time [hours]')
ylabel('Argument of Perigee [degree]')

figure
plot(t, ta)
xlabel('time [hours]')
ylabel('True Anomaly [degree]')



