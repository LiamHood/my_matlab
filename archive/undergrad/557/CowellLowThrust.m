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
        cosu = dot( n , r )/( norm( n )*norm( r ) ) ;
        sinu = norm( cross( n , r ) )/( norm( n )*norm( r ) ) ;
        u = atan2d( sinu , cosu ) ;
        if r(3) < 0
            u = 360 - u ;
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