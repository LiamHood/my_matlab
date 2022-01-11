function [ v1 , v2 ] = Lambert_UV( r1 , r2 , dt , tm , mu )
% uses universal variable to find velocity at two postions 
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
    A = tm*sqrt( r1mag*r2mag*( 1 + cos_ta ) ) ;
    if A >= -1e-12 && A <= 1e-12
        error( 'A is too close to zero' )
    end
    z = 0 ;
    [ c3 , c2 ] = StumpfCalc( z )  ;
    
    zhi = 4*pi^2 ;
    zlow = -4*pi^2 ;
    
    var = 1 ;
    tol = 1e-3 ;
    iteration = 0 ;
    lim = 1e6 ;
    while var > tol
        y = r1mag + r2mag + ( A*( z*c3 - 1 ) )/sqrt( c2 ) ;
        if A > 1e-4 
            if y < -1e-4
                zlow = ( ( ( -r1mag - r2mag )*sqrt(c2)/A ) + 1 )/c3 ;
            end
        end
        chi = sqrt( y/c2 ) ;
        dtn = ( chi^3*c3 + A*sqrt( y ) )/sqrt( mu ) ;
        if dtn <= dt
            zlow = z ;
        else
            zhi = z ;
        end
        z = ( zhi + zlow )/2 ;
        [ c3 , c2 ] = StumpfCalc( z )  ;
        
        var = abs( dtn - dt ) ;
        iteration = iteration + 1 ;
        if iteration > lim
            error( 'dt value failed to converge' )
        end
    end
    
    f = 1 - y/r1mag ;
    gdot = 1 - y/r2mag ;
    g = A*sqrt( y/mu ) ;
        
    v1 = ( r2 - f*r1 )/g ;
    v2 = ( gdot*r2 - r1 )/g ;
end