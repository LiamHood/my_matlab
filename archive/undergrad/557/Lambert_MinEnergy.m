function [ v1 , amin , emin , tmin , tp ] = Lambert_MinEnergy( r1 , r2 , mu )
% finds minimum energy, or fundamental ellipse, between two points
% r1,2 are position vectors, columns, in km
% v1 is column vectors in km/s
% amin is semimajor axis of min energy ellipse
% emin is eccentricity
% tmin is time of transit, short way then long way
% tp is parabolic time of transit

r1mag = norm( r1 ) ;
r2mag = norm( r2 ) ;
cos_ta = dot( r1 , r2 )/( r1mag*r2mag ) ;
c = sqrt( r1mag^2 + r2mag^2 - 2*r1mag*r2mag*cos_ta ) ;
s = ( r1mag + r2mag + c )/2 ;
amin = s/2 ;
pmin = ( r1mag*r2mag / c )*( 1 - cos_ta ) ;
emin = sqrt( 1 - ( 2*pmin )/s ) ;
beta = 2*asin( sqrt( ( s - c )/s ) ) ;
tmin(1) = sqrt( amin^3 / mu )*( pi - ( beta - sin( beta ) ) ) ;
tmin(2) = sqrt( amin^3 / mu )*( pi + ( beta - sin( beta ) ) ) ;
tp = ( 1/3 )*sqrt( 2/mu )*( s^1.5 - ( s - c )^1.5 ) ;
v1 = sqrt( mu*pmin )/( r1mag*r2mag*sin( acos( cos_ta ) ) )*( r2 - ( 1 - ( r2mag/pmin )*( 1 - cos_ta ) )*r1 ) ;
    
end