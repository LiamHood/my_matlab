function [ ap ] = ThreeBody( t , r , JDo , body ) 
rmag = norm( r ) ;
mu = 398600 ;
JD = JDo + t/( 24*60*60 ) ;
n = JD - 2451545.0 ;
if contains( body , "Sun" )
   mu = 132.712e9 ;
   rbavg = 149597870.691 ;
   M = 357.529 + 0.98560023*n ;
   L = 280.459 + 0.98564736*n ;
   lam = L + 1.915*sind( M ) + 0.0200*sind( 2*M ) ;
   eps = 23.439 - ( 3.56e-7 ) * n ;
   us = [ cosd( lam ) ; sind( lam )*cosd( eps ) ; sind( lam )*sin( eps ) ] ;
   rbmag = ( 1.00014 - 0.01671*cosd( M ) - 0.000140*cos( 2*M ) )*rbavg ;
   rb = rbmag*us ;
elseif contains( body , "Moon" )
    [X,Y,Z]=moonpos(JD,'q2000') ; 
    rb = [ X ; Y ; Z ] ;
end
rbmag = norm( rb ) ;
rbs = rb - r ;
q = dot( r , ( 2*rb - r ) )/( rbmag^2 ) ;
Fq = q^2 - 3*q + 3 / ( 1 + ( 1 - q )^1.5 ) * q ;
ap = mu/( norm( rbs )^3 )*( Fq*rb - r ) ;
end