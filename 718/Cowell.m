function [ t,  rv, COES] = Cowell( tspan , r0v0 , mu , tol , forces , A , m, Re)
fprintf( "Starting Cowell \n" )
    opts = odeset('RelTol', tol, 'AbsTol', tol, 'Events', @Reentry) ;
    [ t , rv ] = ode45(@CowellForceFun, tspan, r0v0, opts, mu, forces, A, m, Re) ;
    for ii = 1:length(t)
        COES(ii,:) = state2COE( rv(ii,:) , mu ) ;
    end
    
    function drdv = CowellForceFun(t, rv, mu, forces, A, m, Re)
        t/(3600*24*365)
        r = rv(1:3) ;
        v = rv(4:6) ;
            if contains( forces , "drag" )
                apd = DragAcceleration( r , v , A , m ) ;
            else
                apd = zeros(3,1) ;
            end
            if contains( forces , "gravity" )
                if contains( forces , "J2" )
                    aph = J2accel( r ) ; 
                end
                if contains( forces , "J3" )
                    aph = aph + J3accel( r ) ; 
                end
            else 
                aph = zeros(3,1) ; 
            end
            ap = apd + aph;
        dr = v ;
        dv = ( -mu / norm(r)^3 )*r + ap ;
        drdv = [ dr ; dv ] ;
    end

    function [ value , isterminal , direction ] = Reentry(t, rv, mu, forces, A, m, Re)
        r = norm( rv(1:3) ) ;
        value = r - ( Re + 120 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end
end