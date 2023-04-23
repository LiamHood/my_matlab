clear; close all; clc;

Re = 6378.1;
mu = 398600;

alt = 550;
a0 = alt + Re;
ecc0 = .005 ;
inc0 = sso(a0, ecc0);

s2d = 1/(3600*24*365) ;
r2d = 180/pi ;
d2r = pi/180 ;

raan0 = 0*d2r ;
omega0 = 00*d2r ;
theta0 = 0*d2r ;
p = a0*( 1 - ecc0^2 ) ;
h0 = sqrt( mu*p ) ;

A = .3*.1;
m = 4;

% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
COESo = [h0, inc0, ecc0, raan0, omega0, theta0];
[r0, v0 ] = coes2state(COESo, mu ) ;
r0v0 = [r0; v0 ];

tend = 25/s2d;
tspan = [0, tend];
tol = 1e-8;
dt = 1*3600;
forces = "drag, gravity, J2" ;
[t, state, COES] = Encke(dt, tspan, r0, v0, mu, forces, A, m, tol, Re);
% [t, state, COES] = Cowell(tspan , [r0; v0], mu, tol, forces, A, m, Re);
% [t, state, COES] = VoP(tspan, [r0; v0], mu, tol, forces, A, m, Re);

for ii = 1:length(t)
    r(ii) = norm(state(1:3, ii));
end

figure
axis([t(1), t(end)*s2d, 0, a0*(1 + ecc0) - Re])
hold on
plot(t.*s2d, r-Re, ".")
plot(t.*s2d, COES( :,8 )-Re)
plot(t.*s2d, COES( :,9 )-Re)
plot(t.*s2d, COES( :,7 )-Re)
ylabel( 'Altitude [km]' )
xlabel( 'Time [year]' )
legend( 'Altitude', 'Perigee' , 'Apogee', 'Semi-Major Axis' )
hold off

figure
axis([t(1), t(end)*s2d, 0, 360])
plot(t.*s2d, COES(:,4)./d2r, ".")
ylabel( 'RAAN [degree]' )
xlabel( 'Time [years]' )


% figure
% plot3(state(1,:), state(2,:), state(3,:))

index = max(find(r > Re + 120));
fprintf("The time to deorbit to 120 km is %f years\n", t(index)*s2d)
[T, ecl_time, ecl_proportion] = eclipse(a0, mu, Re);

Je = 237;
Js = 1371;
a = .33;
F = .15;
F_ecl = .7;
A_alpha = .3*.3*pi;
A_eps = .3*.3*pi*4;
alpha = 1;
eps = 1;
Q = 0;
[T_no_ecl, T_ecl] = avg_temp(Re, Je, Js, a0, a, F, F_ecl, ecl_proportion, A_alpha, A_eps, alpha, eps, Q);
