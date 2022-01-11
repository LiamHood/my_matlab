function [ h , inc , ecc , RAAN , omega , theta , a ] = OrbitalElements( r , v , mu )
% all angles output in degrees
    r2d = 180/pi ; % radians to degrees

    Kh = [ 0 0 1 ] ; % K hat

    distance = norm( r ) ;
    speed = norm( v ) ;
    vr = dot( r , v )/distance ; % radial velocity
    h = cross( r , v ) ; % specific angular momentum
    hmag = norm( h ) ; % specific angular momentum
    inc = acos(h(3)/norm(h)) ; %inclination
    eccv = (1/mu)*( cross(v,h)-mu*(r/distance) ) ; %eccentricity vector
    ecc = norm( eccv ) ; % eccentricity
    Nv = cross( Kh , h ) ; % Node line
    N = norm( Nv ) ;

    if Nv(2) > 0 
        RAAN = acos(Nv(1)/N) ; %Right ascension of ascending node
    elseif Nv(2) < 0 
        RAAN = 2*pi - acos(Nv(1)/N) ; %Right ascension of ascending node
    else
        RAAN = 'Undefined' ;
    end

    if eccv(3) > 0 
        omega = acos(dot(Nv,eccv)/(N*ecc)) ; % Argument of perigee
    elseif eccv(3) < 0 
        omega = 2*pi - acos(dot(Nv,eccv)/(N*ecc)) ; % Argument of perigee
    else
        omega = 'Undefined' ;
    end

    % True anomaly
    if vr >= 0
        theta = acos( dot(eccv,r)/(ecc*distance) ) ; 
    else
        theta = 2*pi - acos( dot(eccv,r)/(ecc*distance) ) ; 
    end

    epsilon = speed^2/2 - mu/distance ; % specific energy
    a = - mu/(2*epsilon) ; % semi-major axis
    
    hvec = h ;
    h = hmag ;
    inc = r2d*inc ;
    RAAN = r2d*RAAN ;
    omega = r2d*omega ;
    theta = r2d*theta ;
end


