function HW2_LiamHood()
%% Aero 421 HW2
% Liam Hood
clear; close all; clc;


%% 12.1
% Attitude
roll = pi/4 ;
pitch = pi/4 ;
yaw = pi/4 ;
[ CbG ] = C_321Euler( roll , pitch , yaw ) ;

% position
Ro = 7000e3 ; % m
Ro_G = [ 0 , 0 , Ro ]' ; % ECI
Ro_b = CbG*Ro_G ; % body
Rbcross = crossmatrix( Ro_b ) ;

% Inertia
Ip = zeros( 3,3 ) ; % kg*m^2
Ip(1,1) = 100 ;
Ip(2,2) = 120 ;
Ip(3,3) = 80 ;

mue = 398600 ; % km^3/s^2
mu = 398600*(1000)^3 ; % m^3/s^2

Tgg = (( 3*mu )/( Ro^5 ))*Rbcross*Ip*Ro_b ; 
disp( 'The gravity gradient torque in the body frame in N*m is' )
disp( Tgg )

%% 13.1

I = PrincipalMoI( 60 , 140 , 150 ) ; % kg*m^2
w0 = [ .1 , .02 , .5 ] ; % rad/s
opts = odeset( 'RelTol' , 1e-10 , 'AbsTol' , 1e-10 ) ;
[ t , w ] = ode45( @TorqueFree , [ 0 , 1e3 ] , w0 , opts , I ) ;
figure
plot( t , w(:,1) , t , w(:,2) , t , w(:,3) )
title( 'Angular Velocity' )
xlabel( 'Time (s)' )
ylabel( 'Angular Velocity (rad/s)' )
legend( 'w_x' , 'w_y' , 'w_z' )
disp( 'Angular velocity in z is constant and the others oscillate. Both of ' )
disp( 'the transverse velocities have the same amplitude which is ' )
disp( 'because I didn''t require the s/c to be axisymmetric yet. ' )
disp( ' ' )

disp( 'I assumed axisymmetric for this part because we have only derived the' )
disp( 'equations for nutation and precession rate for that case ' )
hx = I(1,1) .* w(:,1) ;
hy = I(2,2) .* w(:,2) ;
hz = I(3,3) .* w(:,3) ;
ht = ( hx.^2 + hy.^2 ).^(1/2) ;
h = ( hx.^2 + hy.^2 + hz.^2 ).^(1/2) ;

nut = asin(ht./h) ;
nutation = mean( nut ) ;
disp([ 'The nutation angle is ' , num2str(nutation) , ' rad' ])
precessionrate = mean( h./( mean( I(1,1) , I(2,2) )) ) ;
disp([ 'The precession rate is ' , num2str(precessionrate) , ' rad/s' ])
end
