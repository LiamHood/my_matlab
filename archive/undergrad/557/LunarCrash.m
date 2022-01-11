clear ; close all ; clc ;
mu = 4904.8695 ;
Rm = 1737 ;
d2r = pi/180 ;

a = 7500 ;
inc = 55 ;
raan = 0 ;
ecc = 1e-9 ;
omega = 90 ; 
theta = 0 ;
p = a*(1-ecc^2) ;
h = sqrt( mu*p ) ;
COES = [ h , inc*d2r , ecc , raan*d2r , omega*d2r , theta*d2r ] ;
[ r0 , v0 ] = coes2state( COES , mu ) ;
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 

T = 0.014 ; % Newtons
Isp = 1100 ; % s
mi = 150 ; % kg
lat = 0 ; % set
udiff = 30  ;
ucenter = asind( lat/inc ) + 180 ;

    umin = ucenter - udiff ;
    umax = ucenter + udiff ;
    tspan = [ 0 , 2000 ]*86400 ;
    [ tlow , rlow , vlow , mlow ] = CowellLowThrust( tspan , [ r0 ; v0 ; mi ] , mu , 1e-8 , T , Isp , 7000 , 0  , 360 ) ;
    [ tm , rm , vm , mm ] = CowellLowThrust( tspan , [ rlow(:,end) ; vlow(:,end) ; mlow(end) ] , mu , 1e-8 , T , Isp , Rm+50 , umin , umax ) ;
    [ tf , rf , vf , mf ] = CowellLowThrust( tspan , [ rm(:,end) ; vm(:,end) ; mm(end) ] , mu , 1e-8 , T , Isp , Rm , ucenter - 1 , ucenter + 1 ) ;
    t = [ tlow ; tm+tlow(end) ; tf ] ;
    r = [ rlow , rm , rf ] ;
    v = [ vlow , vm , vf ] ;
    m = [ mlow , mm , mf ]' ;
    COESf = state2coes_display( [ r(:,end) ; v(:,end) ] , mu )  ;
    % fprintf( '%f \n' , [ COESf(2:5)' ; COESf(7:9)' ] )
    h = cross( r(:,end) , v(:,end) ) ;
        n = cross( [ 0 ; 0 ; 1 ] , h/norm( h ) ) ;
        cosu = dot( n , r(:,end) )/( norm( n )*norm( r(:,end) ) ) ;
        sinu = norm( cross( n , r(:,end) ) )/( norm( n )*norm( r(:,end) ) ) ;
        uf = atan2d( sinu , cosu ) ;
        if r(3) < 0
            uf = 360 - uf ;
        end
%     mout(ii) = mi - m(end) ;
%     dout(ii) = ( t(end)+tlow(end) )/86400 ;
%     eout(ii) = abs( abs( COESf(2)*sind( uf ) ) - abs( lat ) ) ;
%     eeout(ii) = abs( abs( COESf(2)*sind( COESf(5) ) ) - abs( lat ) ) ;

% figure
% plot( udiff , mout )
% title( 'Mass Used' )
% figure
% plot( udiff , dout )
% title( 'Duration' )
% figure
% plot( udiff , eout )
% title( 'Error in Latitude' )
% figure
% plot( udiff , eeout )
% title( 'Error in Latitude' )
% figure
% plot( t/86400 , m )
figure
axis equal
plot3( r(1,:) , r(2,:) , r(3,:) )

fprintf( 'Range of Orbit Thrust Applied: %f degrees \n' , 2*udiff )
fprintf( 'Latitude of Collision: %f \n' , COESf(2)*sind( uf ) )
fprintf( 'Latitude of Perilune: %f \n' , COESf(2)*sind( COESf(5) ) )
fprintf( 'Planned Latitude: %f \n' , lat )
fprintf( 'Radius of Apolune: %f km \n' , COESf(8) )
fprintf( 'Radius of Perilune: %f km \n' , COESf(9) )
fprintf( 'Duration: %f days \n' , ( t(end)+tlow(end) )/86400 )
fprintf( 'Mass used: %f kg \n' , mi - m(end) )