clear; close all; clc;
Re = 6378.1;
mu = 398600;

alt = 554.283473;
a0 = alt + Re;
ecc0 = .005 ;
inc0 = sso(a0, ecc0);

s2d = 1/(3600*24*365) ;
r2d = 180/pi ;
d2r = pi/180 ;
% lla2eci([40,-102,COES_ref(8)],[2023,6,1,15,12,0])
offsets = eci2lla([Re, 0, 0], [2023,4,25,5,12,0]);
lat0 = 37;
long0 = 95.12658;
additional = (deg2rad(lat0)*Re/tan(inc0))/Re;
GMAT_raan0 = 360-(rad2deg(additional)+long0+offsets(2));
raan0 = 360-(additional+deg2rad(long0));
omega0 = (deg2rad(lat0)*Re/sin(inc0))/Re;
GMAT_omega0 = rad2deg(omega0);
kbps = 1.6e7/413;
theta0 = 0*d2r ;
p = a0*( 1 - ecc0^2 ) ;
h0 = sqrt( mu*p ) ;

A = .3*.1;
m = 4;

% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
COESo = [h0, inc0, ecc0, raan0, omega0, theta0];
[r0, v0 ] = coes2state(COESo, mu ) ;
COES_ref = state2coes([r0;v0], mu);
% lla2eci([40,-102,COES_ref(8)],[2023,6,1,15,12,0])

r0v0 = [r0; v0 ];

tend = 25/s2d;
tspan = [0, tend];
tol = 1e-8;
dt = 1*3600;
forces = "drag" ;
[t, state, COES] = Encke(dt, tspan, r0, v0, mu, forces, A, m, tol, Re);
% [t, state, COES] = Cowell(tspan , [r0; v0], mu, tol, forces, A, m, Re);
% [t, state, COES] = VoP(tspan, [r0; v0], mu, tol, forces, A, m, Re);

for ii = 1:length(t)
    r(ii) = norm(state(1:3, ii));
end

figure
axis([t(1), t(end)*s2d, 0, a0*(1 + ecc0) - Re])
hold on
% plot(t.*s2d, r-Re, ".")
plot(t.*s2d, COES( :,8 )-Re)
plot(t.*s2d, COES( :,9 )-Re)
plot(t.*s2d, COES( :,7 )-Re)
ylabel( 'Altitude [km]' )
xlabel( 'Time [year]' )
legend( 'Perigee' , 'Apogee', 'Semi-Major Axis' )
% 'Altitude',
hold off
% 
figure
axis([t(1), t(end)*s2d, 0, 360])
plot(t.*s2d, COES(:,4)./d2r, ".")
ylabel( 'RAAN [degree]' )
xlabel( 'Time [years]' )


% figure
% plot3(state(1,:), state(2,:), state(3,:))

index = max(find(r > Re + 120));
fprintf("The time to deorbit to 120 km is %f years\n", t(index)*s2d)
fprintf("The inclination is %f degrees\n", rad2deg(inc0))
[T, ecl_time, ecl_proportion] = eclipse(a0, mu, Re);



