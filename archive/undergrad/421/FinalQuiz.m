%% Final Quiz
% Aero 421
% Liam Hood and Jason Dove
clear ; close all ; clc ;

%% Constants
mu = 398600 ; % km^3/s^2 

%% Givens
% COES
h = 53335.2 ; % km^2/s
ecc = 0 ; 
RAAN = 0 ;
inc = 98.43 ; % degrees
omega = 0 ;
theta = 0 ; 

[ r , v ] = coes2state(  h , ecc , theta , RAAN , omega , inc , mu ) ;

angvelb = [ 0 -0.001047 0 ]' ;
% body to LVLH
eps = [ 0 0 0 ] ;
eta = 1 ;
quat0b = [ eta , eps ]' ;
% LVLH to ECI
rd = r/norm(r) ;
vd = v/norm(v) ;
zLVLH = -rd ;
yLVLH = -cross( rd , vd ) ;
xLVLH = cross( yLVLH , zLVLH ) ;
Flvlh = [ xLVLH , yLVLH , zLVLH ] ;
x = [ 1 0 0 ]' ;
y = [ 0 1 0 ]' ;
z = [ 0 0 1 ]' ;
Feci = [ x , y , z ] ;
Cbg0 = Flvlh/Feci ;
quat0 = dcm2quat( Cbg0 )' ;
[ xr , yr , zr ] = dcm2angle( Cbg0 ) ;

% specs
side = 2 ; % side length in m
busmass = 500 ; % mass of bus in kg
magdiz = [ 0 0 -0.5 ] ; % magnetic dipole moment A*m^2

%% Moments of Inertia

% Moments of inertia of bus
Ixbus = (busmass*side^2)/6 ;
Iybus = (busmass*side^2)/6 ;
Izbus = (busmass*side^2)/6 ;

sensmass = 100 ; % sensor mass
% moments of inertia of bus
Izsens = (sensmass*.25^2)/6 ;
Iysens = (1/12)*sensmass*(.25^2+1) + sensmass*(.5+.25)^2 ;
Ixsens = (1/12)*sensmass*(.25^2+1) + sensmass*(.5+.25)^2 ;

panmass = 20 ; % solar panel mass
% moments of inertia of panel
Iypan = 2*(1/12)*panmass*(.05^2+2^2) ;
Izpan = 2*( (1/12)*panmass*(2^2+3^2) + panmass*(2^2) ) ;
Ixpan = 2*( (1/12)*panmass*(.05^2+3^2) + panmass*(2^2) ) ;

% moments of inertia of s/c
Ix = Ixbus + Ixsens + Ixpan ;
Iy = Iybus + Iysens + Iypan ;
Iz = Izbus + Izsens + Izpan ;
Ip = [ Ix 0 0 ; 0 Iy 0 ; 0 0 Iz ] ;

% Areas
Ax = .05*3*2 + 2*2 + .25*1 ;
Ay = 2*2 + .25*1 ;
Az = 3*2*2 + 2*2 ;
Avec = [ Ax , Ay , Az ] ;

%% Run Simulation
opts = odeset( 'Reltol' , 1e-8 , 'Abstol' , 1e-8 ) ;
tspan = [ 0 , 11*60*60 ] ;
UTC = [ 2019 , 3 , 19 , 12 , 0 , 0 ] ;
JDay = juliandate( UTC ) ;
j2000 = juliandate( [ 2000 , 1 , 1 ] ) ;
Day0 = JDay - j2000 ;

state0 = [ quat0 ; angvelb ; r ; v ] ; 
[ t , state ] = ode45( @Simulation , tspan , state0 , opts , Ip , Avec , Day0 , magdiz , UTC ) ;

figure
plot( t , state(:,5) , t , state(:,6) , t , state(:,7) )
title( 'Angular Velocity' )
xlabel( 'Time' )
ylabel( 'Angular Velocity (rad/s)' )
legend( 'x' , 'y' , 'z' )

for ii = 1:length( t )
    quat(ii,:) = state( ii , 1:4 ) ;
    Cbg = quat2dcm( quat(ii,:) ) ;
    [ pitch(ii) , roll(ii) , yaw(ii) ] = dcm2angle( Cbg ) ;
    quatadj(ii,:) = dcm2quat( Cbg ) ;
end

figure
plot( t , pitch*(180/pi) , t , roll*(180/pi) , t , yaw*(180/pi) )
title( 'Angular Position' )
xlabel( 'Time' )
ylabel( 'Angle (degrees)' )
legend( 'Roll' , 'Pitch' , 'Yaw' )

figure
plot( t , quatadj(:,1) , t , quatadj(:,2) , t , quatadj(:,3) , t , quatadj(:,4) )
title( 'Quaternions over Time' )
xlabel( 'Time' )
ylabel( 'Quaternion Value' )
legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )