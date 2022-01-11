%% Homework 1
% Liam Hood
% Aero 351

function Aero_351_HW1

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
disp( ' ' )

%% Part 2 

disp( '2' )

% a
disp( 'a' )
    time_a = [ 10 , 0 , 0 ] ; %Time in hours,minutes,seconds
    date_a = [ 21 , 12 , 2007 ] ; %Date in day,month,year
    loc_a = 144+58/60 ; %Location in east longitude
    [ lst_a ] = LSidereal( time_a , date_a , loc_a) ;
    disp([ 'The LST at ' , num2str(loc_a) , ' degrees east at ' , num2str(date_a) , ' at ' , num2str(time_a(1)) , ':' , num2str(time_a(2)) , ':' , num2str(time_a(3)) , ' is '])
    disp([ num2str(lst_a ) , ' degrees' ]) %display answer
% b
disp( 'b' )
    time_b = time ; %Time in hours,minutes,seconds
    date_b = [ 4 , 7 , 2018 ] ; %Date in day,month,year
    loc_b = -120.653 ; %location in east longitude
    [ lst_b ] = LSidereal( time_b , date_b , loc_b) ;
    disp([ 'The LST at ' , num2str(loc_b) , ' degrees east at ' , num2str(date_b) , ' at ' , num2str(time_b(1)) , ':' , num2str(time_b(2)) , ':' , num2str(time_b(3)) , ' is ' ])
    disp([ num2str(lst_b) , ' degrees' ]) %Display answer
    disp( ' ' )
    
%% Part 3

disp( '3. Curtis 2.4' )

mu_e = 398600 ; % km^3/s
rv = [ 3207 5459 2714 ] ; %initial position (km)
vv = [ -6.532 .7835 6.142 ] ; %initial velocity (km/s)
m2 = 1000 ; %mass sattelite
t0 = 0 ; %initial time
tf = 24*60*60 ; %final time in seconds

istate = [ rv vv ] ; % state vectors

ooptions = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[ ~ , nstate ] = ode45( @TwoBodyMotion , [ t0 tf ] , istate , ooptions , mu_e ) ;

figure
plot3( nstate(:,1) , nstate(:,2) , nstate(:,3) )
xlabel( 'x (km)' )
ylabel( 'y (km)' )
zlabel( 'z (km)' )

disp([ 'Position after 24 hours is ' , num2str([ nstate(1) nstate(2) nstate(3) ]) , 'km' ])
disp([ 'Radius is ' , num2str( norm([ nstate(1) nstate(2) nstate(3) ]) ) , ' km' ])
disp([ 'Velocity after 24 hours is ' , num2str([ nstate(4) nstate(5) nstate(6) ]) , ' km/s' ])
disp([ 'Speed is ' , num2str( norm([ nstate(4) nstate(5) nstate(6) ]) ) , ' km/s' ])
disp( ' ' )

%% Part 4

disp( '4. Curtis 2.10' )
disp( 'Start with velocity component equations' )
disp( 'v=sqrt(v_r^2+v_az^2) and v_az=(mu/h)*(1+e*cos(theta)) and v_r=(mu/h)*e*sin(theta)' )
disp( 'Then combine into single v equation and simplify' )
disp( 'v=sqrt( ((mu/h)*(1+e*cos(theta)))^2 + ((mu/h)*e*sin(theta))) ' )
disp( 'v=sqrt( (mu/h)^2*(1+2*e*cos(theta)+e^2*cos(theta)^2+e^2*sin(theta)^2))' )
disp( 'v=(mu/h)*sqrt(1+2*e*cos(theta)+e^2*(cos(theta)^2+sin(theta)^2 ' )
disp( 'v=(mu/h)*sqrt(1+2*e*cos(theta)+e^2)' )
disp( ' ' )


%% Part 5

disp( '5. Curtis 2.17' )

mu_m = 42828.37 ; % mu of mars (km^s/s)
r_m = 3390 ; %km
altitude = 200 ; %km
r = r_m + altitude ; %km

period = ( ( 2*pi )/( sqrt(mu_m) ) )*r^1.5 ; % period of orbit in seconds
period_hr = period/(60*60) ; %Period in hours
speed = sqrt( mu_m/r ) ; % Speed of orbit 

disp([ 'The period of the orbit is ' , num2str( period_hr ) , ' hours' ])
disp([ 'The speed of the orbit is ' , num2str( speed ) , ' km/s' ])
disp( ' ' )

%% Part 6 
r_p = 10000 ; % Radius at perigee (km)
r_a = 100000 ; % Radius at apogee (km)
r_e = 6378 ; % Radius of Earth (km)

disp( '6. Curtis 2.21' )

% a 
disp( 'a' )
ecc6 = ( r_a - r_p )/( r_a + r_p ) ; %eccentricity of orbit
disp([ 'The eccentricity is ' , num2str( ecc6 ) ])

h = sqrt( r_p * mu_e * ( 1 + ecc6 ) ) ; %angular momentum

disp( ' ' )
% b
disp( 'b' )
a = ( r_a + r_p )/2 ; % semi-major axis (km)
disp([ 'The semi-major axis is ' , num2str( a ) , ' km' ])
disp( ' ' )

% c 
disp( 'c' )
T = ((( 2*pi )/sqrt(mu_e)) * a^1.5 )/60^2; %Period in hours
disp([ 'The period is ' , num2str( T ) , ' hours' ])
disp( ' ' )

% d
disp( 'd' )
se = - mu_e/( 2*a ) ; %Specific energy of orbit
disp([ 'The specific energy of the orbit is ' , num2str( se ) , ' km^2/s^2' ])
disp( ' ' )

% e 
disp( 'e' )
r6e = 10000 + r_e ; % Radius at altitude of 10000 km
ta6e = acosd( (( h^2/( mu_e*r6e )) - 1 ) / ecc6 ) ; %True anomaly above radius (degrees)
disp([ 'The true anomaly when the satellite is at an altitude of 10,000 km is ' ])
disp([ num2str( ta6e ) , ' degrees' ])
disp( ' ' )

% f
disp( 'f' )
v_az = h/r6e ; % Azmuthal velocity (km/s)
v_r = ( mu_e/h ) * ecc6 * sin( ta6e ) ; %Radial velocity (km/s)
disp([ 'The azmuthal velocity is ' , num2str( v_az ) , ' km/s' ])
disp([ 'The radial velocity is ' , num2str( v_r ) , ' km/s' ]) 
disp( ' ' )

% g
disp( 'g' )
v_per = h/r_p ; % Velocity at perigee (km/s)
v_apo = h/r_a ; % Velocity at apogee (km/s)
disp([ 'The velocity at perigee ' , num2str( v_per ) ])
disp([ 'The velocity at apogee ' , num2str( v_apo ) ])
disp( ' ' )

%% Functions

function [ Jd , Jo , UT] = Julian( time , date )
%Calculates the Julian Date from a date and time
%   Uses an input of date in form [dd,mm,yyyy] and time in UT [hour,minute,second] to find Julian date. BCE years should be
%   negative

% Julian date without time
    Jo = 367*date(3) - floor(( 7*( date(3)+floor(( date(2)+9 )/12 )) )/4)+floor((275*date(2))/9) + date(1) + 1721013.5 ;
% Time 
    hour = time(1) ; %hours past noon as fraction of a day
    minute = time(2)/(60) ; %minutes as fraction of a day
    second = time(3)/(60*60) ; %seconds as fraction of a day
    UT = hour + minute + second ; % add time together 
    Jd = Jo + UT/24 ; % Full Julian date
    
end

function [Lst] = LSidereal( time , date , location )
%Calculates Local Siderial Time 
%   Input arguments are time in [hour,minute,second] and date in
%   [day,month,year] and location is in east longitude. West longitude is
%   negative
J2000 = 2451545 ;
JC = 36525 ;
    [ ~ , Jo , UT ] = Julian( time , date ) ; %Runs Julian to find Jo
    To = ( Jo - J2000 )/JC ;
    theta_G0 = 100.4606184 + 36000.77004*To + .000387933*To^2-2.58e3*To^3 ;
    if theta_G0 > 360 % If to big to 
        while theta_G0 > 360
            theta_G0 = theta_G0 - 360 ;
        end
    elseif theta_G0 < 0
        while theta_G0 < 0
            theta_G0 = theta_G0 + 360 ;
        end
    end
    theta_G = theta_G0 + 360.98564724*(UT/24) ;
    theta = theta_G + location ;
        if theta > 360
        while theta > 360
            theta = theta - 360 ;
        end
    elseif theta < 0
        while theta < 0
            theta = theta + 360 ;
        end
        end
        Lst = theta ;
    
end

function dstate_dt = TwoBodyMotion( t , state , mu ) 
% Finds change in state with respect to time. Input time, t, in seconds and
% state as position vector followed by velocity vector as well as mu

rad = norm( [ state(1) state(2) state(3) ] ) ; %radius

dx = state(4) ; % velocity in x
dy = state(5) ; % velocity in y
dz = state(6) ; % velocity in z
ddx = -mu*state(1)/rad^3 ; % acceleration in x
ddy = -mu*state(2)/rad^3 ; % acceleration in y
ddz = -mu*state(3)/rad^3 ; % acceleration in z

dstate_dt = [ dx ; dy ; dz ; ddx ; ddy ; ddz ] ;

end

end
