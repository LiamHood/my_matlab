function [ dstate ] = Simulation( t , state , Ip , Avec , Day0 , mbar , UTC )
% Simulates an attitude of s/c in LEO. State is quaternion [ eta , eps ]
% followed by angular velocity

% to track progress
t 
% Set up
    mu = 398600 ;
    dstate = zeros(13,1) ;
% Set up names
    Ix = Ip(1,1) ;
    Iy = Ip(2,2) ;
    Iz = Ip(3,3) ;
    Ax = Avec(1) ;
    Ay = Avec(2) ;
    Az = Avec(3) ;

    quat = state( 1:4 )' ;
    w = state( 5:7 ) ;
    r = state( 8:10 ) ;
    v = state( 11:13 ) ;
    
    eta = quat( 1 ) ;
    eps = quat( 2:4 ) ;
    Cbg = quat2dcm( quat ) ;
    
    
     % FIX THIS
     if t >= 60
         tmin = floor( t/60 ) ;
         t = t - tmin*60 ;
         if tmin >= 60
             thour = floor( tmin/60 ) ;
             tmin = tmin - thour*60 ;
         else
             thour = 0 ;
         end
     else
         tmin = 0 ;
         thour = 0 ;
     end
     if ( thour + UTC(4) ) >= 24 
         thour = thour - 24 ;
         UTC(3) = UTC(3) + 1 ;
     end
     UTCnow = UTC + [ 0 0 0 thour tmin t ] ;

    JDay = Day0 + t/(24*60*60) ; % Find Julian day of current moment

% Find change in quaternions
    epscross = crossmatrix( eps ) ;
    epsdot = .5 * ( eta*eye(3) + epscross ) * w ;
    etadot = .5 * eps * w ;
    dstate(1:4) = [ etadot , epsdot' ]' ;
    

s = [ 1 0 0 ] ;
% Find torques
    Tsrp = SolarRadiationPressure( Cbg , s , Ax , Ay , Az ) ;
    Tgg = GravityGradientTorque( Cbg , Ix , Iy , Iz , r ) ;
    Tmag = MagneticTorque( Cbg , mbar , r , UTCnow ) ;
    Tatmo = AtmoTorque( Cbg , v , Ax , Ay , Az ) ;
    
    T = Tsrp + Tgg + Tmag + Tatmo ;
    
% Find change in angular velocity from torques 
    wdot = zeros( 3,1 ) ;
    wdot(1) = (T(1)-( Iz - Iy )*w(2)*w(3))/Ix ;
    wdot(2) = (T(2)-( Ix - Iz )*w(1)*w(3))/Iy ;
    wdot(3) = (T(3)-( Iy - Ix )*w(2)*w(1))/Iz ;
    dstate(5:7) = wdot ;
    
% New ECI state vectors
    rad = norm( [ r(1) r(2) r(3) ] ) ; %radius

    dstate(8) = v(1) ; % velocity in x
    dstate(9) = v(2) ; % velocity in y
    dstate(10) = v(3) ; % velocity in z
    dstate(11) = -mu*r(1)/rad^3 ; % acceleration in x
    dstate(12) = -mu*r(2)/rad^3 ; % acceleration in y
    dstate(13) = -mu*r(3)/rad^3 ; % acceleration in z


    