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
udiff = 10:2.5:180  ;
ucenter = asind( lat/inc ) + 180 ;

for ii = 1:length( udiff )
    umin = ucenter - udiff(ii) ;
    umax = ucenter + udiff(ii) ;
    tspan = [ 0 , 5000 ]*86400 ;
    [ tlow , rlow , vlow , mlow ] = CowellLowThrust( tspan , [ r0 ; v0 ; mi ] , mu , 1e-10 , T , Isp , 7000 , 0  , 360 ) ;
    [ tm , rm , vm , mm ] = CowellLowThrust( tspan , [ rlow(:,end) ; vlow(:,end) ; mlow(end) ] , mu , 1e-10 , T , Isp , Rm+100 , umin , umax ) ;
    [ tf , rf , vf , mf ] = CowellLowThrust( tspan , [ rm(:,end) ; vm(:,end) ; mm(end) ] , mu , 1e-10 , T , Isp , Rm , ucenter - 1 , ucenter + 1) ;
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
    mout(ii) = mi - m(end) ;
    dout(ii) = ( t(end)+tlow(end) )/86400 ;
    eout(ii) = abs( abs( COESf(2)*sind( uf ) ) - abs( lat ) ) ;
    eeout(ii) = abs( abs( COESf(2)*sind( COESf(5) ) ) - abs( lat ) ) ;
end
figure
plot( udiff*2 , mout )
title( 'Mass Used' )
xlabel( 'Degrees of Orbital Arc Thrust occurs over' )
ylabel( 'Mass Used [kg]' )
figure
plot( udiff*2 , dout )
title( 'Duration' )
xlabel( 'Degrees of Orbital Arc Thrust occurs over' )
ylabel( 'Time to Deorbit [days]' )
figure
plot( udiff*2 , eout )
title( 'Error in Latitude from Target' )
xlabel( 'Degrees of Orbital Arc Thrust occurs over' )
ylabel( 'Latitude Error' )
figure
plot( udiff*2 , eeout )
title( 'Error in Latitude from Perilune' )
xlabel( 'Degrees of Orbital Arc Thrust occurs over' )
ylabel( 'Latitude Error' )
% figure
% plot( t/86400 , m )
% figure
% plot3( r(1,:) , r(2,:) , r(3,:) )
% fprintf( 'Range of Orbit Thrust Applied: %f degrees \n' , 2*udiff )
% fprintf( 'Latitude of Collision: %f \n' , COESf(2)*sind( uf ) )
% fprintf( 'Latitude of Perilune: %f \n' , COESf(2)*sind( COESf(5) ) )
% fprintf( 'Planned Latitude: %f \n' , lat )
% fprintf( 'Radius of Apolune: %f km \n' , COESf(8) )
% fprintf( 'Radius of Perilune: %f km \n' , COESf(9) )
% fprintf( 'Duration: %f days \n' , ( t(end)+tlow(end) )/86400 )
% fprintf( 'Mass used: %f kg \n' , mi - m(end) )

function [ t , r , v , m ] = CowellLowThrust( tspan , r0v0m , mu , tol , T , Isp , Rcb , umin , umax )
%[ t , r , v ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    [ t , rv ] = ode45( @CowellForceFun , tspan , r0v0m , opts , mu , T , Isp , Rcb , umin , umax ) ;
    r = rv(:,1:3)' ;
    v = rv(:,4:6)' ;
    m = rv(:,7)' ;
    
    function drdv = CowellForceFun( t , rvm , mu , T , Isp , Rcb , umin , umax )
        r = rvm(1:3) ;
        v = rvm(4:6) ;
        m = rvm(7) ;
        h = cross( r , v ) ;
        n = cross( [ 0 ; 0 ; 1 ] , h/norm( h ) ) ;
        u = acosd( dot( n , r )/( norm( n )*norm( r ) ) ) ;
        if r(3) < 0 
            u = u + 180 ;
        end
        dr = v ;
        if u >= umin && u <= umax
            aT = ( T/m )*1e-3 ;
            dm = -T/(Isp*9.81) ;
        else
            aT = 0 ;
            dm = 0 ;
        end
        dvu = ( -mu / norm(r)^3 )*r ;
        dv = dvu - aT*v/norm(v) ;
        drdv = [ dr ; dv ; dm ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , rv , mu , T , Isp , Rcb , umin , umax )
        r = norm( rv(1:3) ) ;
        value = r - Rcb ;
        isterminal = 1 ;
        direction = 0 ;
    end
end

function [ r , v ] = coes2state( COES , mu )
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
h = COES(1) ;
inc = COES(2) ;
ecc = COES(3) ;
RAAN = COES(4) ;
omega = COES(5) ;
theta = COES(6) ;

    r_peri = (h^2/mu) * ( 1/( 1 + ecc*cos(theta) ) ) * [ cos( theta ) ; sin( theta ) ; 0 ] ;
    v_peri = (mu/h) * [ -sin( theta ) ; ecc+cos(theta) ; 0 ] ;

    Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
    Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
    Q(1,3) = sin(RAAN)*sin(inc) ;
    Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
    Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
    Q(2,3) = -cos(RAAN)*sin(inc) ;
    Q(3,1) = sin(inc)*sin(omega) ;
    Q(3,2) = sin(inc)*cos(omega) ;
    Q(3,3) = cos(inc) ;

    r = Q*r_peri ;
    v = Q*v_peri ;
end

function COES = state2coes_display( state, mu )
% all angles output in degrees
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
r2d = 180/pi ;
    
    Kh = [ 0 0 1 ] ; % K hat
    R = state(1:3) ;
    V = state(4:6) ;
    r = norm( R ) ;
    v = norm( V ) ;
    vr = dot( R , V )/r ; % radial velocity
    H = cross( R , V ) ; % specific angular momentum
    h = norm( H ) ; % specific angular momentum
    inc = acos( H(3)/h ) ; %inclination
    ECC = (1/mu)*( ( v^2 - mu/r )*R - r*vr*V ) ; %eccentricity vector
    ecc = norm( ECC ) ; % eccentricity
    N = cross( Kh , H ) ; % Node line
    n = norm( N ) ;
if n ~= 0
    RAAN = acos(N(1)/n) ; %Right ascension of ascending node
    if N(2) < 0 
        RAAN = 2*pi - RAAN ; %Right ascension of ascending node
    end
else
    RAAN = 0 ;
end 
    
if n ~= 0 
    if ecc >= 0 
        omega = acos(dot(N,ECC)/(n*ecc)) ; % Argument of perigee
        if ECC(3) < 0 
            omega = 2*pi - omega ; % Argument of perigee
        end
    else
        omega = 0 ;
    end
else
    omega = 0 ;
end
    
if ecc > 0
    theta = acos( dot( ECC , R )/( ecc*r ) ) ;         
    if vr < 0 
        theta = 2*pi - theta ; 
    end
else
    cp = cross( N , R ) ;
    if cp(3) >= 0
        theta = acos( dot( N , R )/( n*r ) ) ;
    else
        theta = 2*pi - acos( dot( N , R )/( n*r ) ) ;
    end
end

    a = (h^2)/( mu*( 1 - ecc^2 ) ) ; % semi-major axis
    rp = a*( 1 - ecc ) ;
    ra = a*( 1 + ecc ) ;

    
    COES = [ h , inc*r2d , ecc , RAAN*r2d , omega*r2d , theta*r2d , a , rp , ra ] ;
end


