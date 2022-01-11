function om_constants_moon

% astrodynamic and utility constants

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dtr rtd mu req j2 egrav

% angular conversion factors

dtr = pi / 180.0;

rtd = 180.0 / pi;

% earth gravitational constant (km**3/sec**2)

mu = 4902.799;

% earth equatorial radius (kilometers)

req = 1738;

% earth oblateness gravity coefficient (non-dimensional)

j2 = 0.0002027 ;

% earth surface gravity (meters/second**2)

egrav = 1.623;


