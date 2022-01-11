clear ;
clc ;
close all ;

r = [ -2429.1 ; 4555.1 ; 4577 ] ; % km in ECI
v = [ -4.7689 ; -5.6113 ; 3.0535 ] ; % km/s in ECI
[ h , inc , ecc , RAAN , omega , theta , a ] = OrbitalElements( r , v , 398600 ) ; % Find orbital elements
