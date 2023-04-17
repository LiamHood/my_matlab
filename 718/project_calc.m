clear; close all; clc;

Re = 6378.1;
mu = 398600;

s2d = 1/(3600*24*365) ;
r2d = 180/pi ;
d2r = pi/180 ;

raan0 = 30*d2r ;
inc0 = 30*d2r ;
omega0 = 50*d2r ;
theta0 = 0*d2r ;
alt = 500;
a0 = alt + Re;
ecc0 = .01 ;
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
dt = 2*3600;
forces = "drag" ;
[t, state, COES] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m, tol, Re);

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
xlabel( 'Time [days]' )
legend( 'Altitude', 'Perigee' , 'Apogee', 'Semi-Major Axis' )
hold off

figure
plot3(state(1,:), state(2,:), state(3,:))

index = max(find(r > Re + 100));
fprintf("The time to deorbit to 150 km is %f years\n", t(index)*s2d)
[T, ecl_time, ecl_proportion] = eclipse(a0, mu, Re);
