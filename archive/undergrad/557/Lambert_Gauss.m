function [ v1 , v2 ] = Lambert_Gauss( r1 , r2 , dt , tm , mu )
% uses gauss method to solve lamberts problem
% best for closely spaced positions
% r1,2 are position vectors, columns, in km
% dt is time of flight in seconds
% tm is short or long way, 1 for short or -1 for long
% v1,v2 are column vectors in km/s
    r1mag = norm( r1 ) ;
    r2mag = norm( r2 ) ;
    cos_ta = dot( r1 , r2 )/( r1mag*r2mag ) ;
    sin_ta = tm*sqrt( 1 - cos_ta^2 ) ;
    if sin_ta > 0
        ta = acos( cos_ta ) ;
    else
        ta = 2*pi + asin( sin_ta ) ;
    end
    cos_hta = cos( ta/2 ) ;
    L = ( r1mag + r2mag )/( 4*sqrt( r1mag*r2mag )*cos_hta ) - .5 ;
    m = ( mu*dt^2 )/( ( 2*sqrt( r1mag*r2mag )*cos_hta )^3 ) ;

    yv = [ 1 , 1 ] ;
    var = 1 ;
    tol = 1e-6 ;
    iteration = 0 ;
    lim = 1e3 ;
    while var >= tol
        yv(1) = yv(2) ;
        x1 = ( m / ( yv(2)^2 ) ) - L ;
%         dE = 2*acos( 1 - 2*x1 ) ;
%         x2 = ( dE - sin( dE ) )/( sin( dE/2 )^3 ) ;
        x2 = ( 4/3 )*( 1 + ( (6*x1)/(5) ) + ( (6*8*x1^2)/(5*7) ) + ( (6*8*10*x1^3)/(5*7*9) ) + ((6*8*10*12*x1^4)/(5*7*9*11))  + ((6*8*10*12*14*x1^5)/(5*7*9*11*13)) ) ;
        yv(2) = 1 + x2*( L + x1 ) ;
        var = abs( yv(2) - yv(1) ) ;
        iteration = iteration + 1 ;
        if iteration > lim
            error( 'y value failed to converge' )
        end
    end
    y = yv(2) ;
    x1 = m/y^2 - L ;
    dE = 2*acos( 1 - 2*x1 ) ;
    a = ( ( dt*sqrt( mu ) )/( 2*y*sqrt( r1mag*r2mag )*cos( ta/2 )*sin( dE/2 ) ) )^2 ;    
    p = ( r1mag*r2mag*( 1 - cos_ta ) )/( r1mag + r2mag - 2*sqrt( r1mag*r2mag )*cos_hta*cos( dE/2 ) ) ;
    % p = ( y^2*r1mag^2*r2mag^2*sin_ta^2 )/( mu*dt^2 ) ;
    f = 1 - ( r2mag/p )*( 1 - cos_ta ) ;
    g = r2mag*r1mag*sin_ta/sqrt( mu/p ) ;
    gdot = 1 - ( r1mag/p )*( 1 - cos_ta ) ;
    % a = ( ( dt*sqrt( mu ) )/( 2*y*sqrt( r1mag*r2mag )*cos_hta*sin( dE ) ) )^2 ;
    f = 1 - ( a/r1mag )*( 1 - cos( dE ) ) ;
    g = dt - sqrt( a^3 / mu )*( dE - sin( dE ) ) ;
    gdot = 1 - ( a/r2mag )*( 1 - cos( dE ) ) ;

    v1 = ( r2 - f*r1 )/g ;
    v2 = ( gdot*r2 - r1 )/g ;
    
end