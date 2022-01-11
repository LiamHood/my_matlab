%% Homework 1
% Liam Hood
% Aero 351

%% Clean up
clear ;
clc ;
close all ;

%% Part 1

disp( '1' )

time = [ 12 , 0 , 0 ] ; %Time in hours,minutes,seconds
date = [ 22 , 9 , 2018 ] ; %Date in day,month,year
[ Jd ] = Julian( time , date ) ; %Finds Julian date
disp([ 'Julian Date of ' , num2str(date) , ' at ' , num2str(time(1)) , ':' , num2str(time(2)) , ':' , num2str(time(3)) , ' UT is' ])
disp( Jd ) %display answer

%% Part 2 

disp( '2' )

% a
disp( 'a' )
    time_a = [ 10 , 0 , 0 ] ; %Time in hours,minutes,seconds
    date_a = [ 21 , 12 , 2007 ] ; %Date in day,month,year
    loc_a = 144+58/60 ; %Location in east longitude
    [ lst_a ] = LSidereal( time_a , date_a , loc_a) ;
    disp([ 'The LST at ' , num2str(loc_a) , ' degrees east at ' , num2str(date_a) , ' at ' , num2str(time_a(1)) , ':' , num2str(time_a(2)) , ':' , num2str(time_a(3)) , ' is '])
    disp( lst_a ) %display answer
% b
disp( 'b' )
    time_b = time ; %Time in hours,minutes,seconds
    date_b = [ 4 , 7 , 2018 ] ; %Date in day,month,year
    loc_b = -120.653 ; %location in east longitude
    [ lst_b ] = LSidereal( time_b , date_b , loc_b) ;
    disp([ 'The LST at ' , num2str(loc_b) , ' degrees east at ' , num2str(date_b) , ' at ' , num2str(time_b(1)) , ':' , num2str(time_b(2)) , ':' , num2str(time_b(3)) , ' is ' ])
    disp( lst_b ) %Display answer
    
%% Part 3

disp( '2.4' )

mu_e = 398600 ; % km^3/s
rv = [ 3207 5459 2714 ] ; %initial position (km)
vv = [ -6.532 .7835 6.142 ] ; %initial velocity (km/s)
m2 = 1000 ; %mass sattelite
t0 = 0 ; %initial time
tf = 24*60*60 ; %final time in seconds

istate = [ rv vv ] ; % state vectors

ooptions = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[ nt , nstate ] = ode45( @TwoBodyMotion , [ t0 tf ] , istate , ooptions , mu_e ) ;

figure
plot3( nstate(:,1) , nstate(:,2) , nstate(:,3) )
xlabel( 'x (km)' )
ylabel( 'y (km)' )
zlabel( 'z (km)' )

disp([ 'Position after 24 hours is ' , num2str([ nstate(1) nstate(2) nstate(3) ]) , 'km' ])
disp([ 'Radius is ' , num2str( norm([ nstate(1) nstate(2) nstate(3) ]) ) , ' km' ])
disp([ 'Velocity after 24 hours is ' , num2str([ nstate(4) nstate(5) nstate(6) ]) , ' km/s' ])
disp([ 'Speed is ' , num2str( norm([ nstate(4) nstate(5) nstate(6) ]) ) , ' km/s' ])

%% Part 5

disp( '2.17' )

mu_m = 42828.37 ; % mu of mars (km^s/s)
r_m = 3390 ; %km
altitude = 200 ; %km
r = r_m + altitude ; %km

period = ( ( 2*pi )/( sqrt(mu_m) ) )*r^1.5 ; % period of orbit in seconds
speed = sqrt( mu_m/r ) ; % Speed of orbit 

disp([ 'The period of the orbit is ' , num2str( period ) , ' km' ])
disp([ 'The speed of the orbit is ' , num2str( speed ) , ' km/s' ])

%% Part 6 
r_p = 10000 ; % Radius at perigee (km)
r_a = 100000 ; % Radius at apogee (km)
r_e = 6378 ; % Radius of Earth (km)

disp( '2.21' )

% a 
disp( 'a' )
ecc6 = ( r_a - r_p )/( r_a + r_p ) ; %eccentricity of orbit
disp([ 'The eccentricity is ' , num2str( ecc6 ) ])

h = sqrt( r_p * mu_e * ( 1 + ecc6 ) ) ; %angular momentum

% b
disp( 'b' )
a = ( r_a + r_p )/2 ; % semi-major axis (km)
disp([ 'The semi-major axis is ' , num2str( a ) , ' km' ])

% c 
disp( 'c' )
T = ((( 2*pi )/sqrt(mu_e)) * a^1.5 )/60^2; %Period in hours
disp([ 'The period is ' , num2str( T ) , ' hours' ])

% d
disp( 'd' )
se = - mu_e/( 2*a ) ; %Specific energy of orbit
disp([ 'The specific energy of the orbit is ' , num2str( se ) , ' km^2/s^2' ])

% e 
disp( 'e' )
r6e = 10000 + r_e ; % Radius at altitude of 10000 km
ta6e = acosd( (( h^2/( mu_e*r6e )) - 1 ) / ecc6 ) ; %True anomaly above radius (degrees)
disp([ 'The true anomaly when the satellite is at an altitude of 10,000 km is ' ])
disp([ num2str( ta6e ) , ' degrees' ])

% f
disp( 'f' )
v_az = h/r6e ; % Azmuthal velocity (km/s)
v_r = ( mu_e/h ) * ecc6 * sin( ta6e ) ; %Radial velocity (km/s)
disp([ 'The azmuthal velocity is ' , num2str( v_az ) , ' km/s' ])
disp([ 'The radial velocity is ' , num2str( v_r ) , ' km/s' ]) 

% g
disp( 'g' )
v_per = h/r_p ; % Velocity at perigee (km/s)
v_apo = h/r_a ; % Velocity at apogee (km/s)
disp([ 'The velocity at perigee ' , num2str( v_per ) ])
disp([ 'The velocity at apogee ' , num2str( v_apo ) ])
