function [ t , r , v] = CowellJ2J3( tspan , r0 , v0 , mu )
%[ t , r , v ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m )
tol = 1e-10 ;
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    
    [ t , rv ] = ode45( @CowellForceFun , tspan , [ r0 ,v0 ] , opts , mu ) ;
    r = rv( : , 1:3 ) ;
    v = rv( : , 4:6 ) ;
    
    function drdv = CowellForceFun( t , rv , mu  )

        r = rv(1:3) ;
        v = rv(4:6) ;
            ap = J2accel( r ) + J2accel( r ) ;
        dr = v ;
        dv = ( -mu / norm(r)^3 )*r + ap ;
        drdv = [ dr ; dv ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , rv , mu )
        r = norm( rv(1:3) ) ;
        value = r - ( 6378 + 100 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end

        function aP = J2accel( r )
            mu = 398600 ;
            J2 = -1.08262617385222e-3 ;
            Re = 6378 ;
            aP = zeros(3,1) ;
            rmag = norm(r) ;
            aP(1) = -3*J2*mu*Re^2*r(1)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
            aP(2) = -3*J2*mu*Re^2*r(2)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
            aP(3) = -3*J2*mu*Re^2*r(3)/(2*rmag^5)*(3-5*r(3)^2/rmag^2) ;
        end
    
        function aP = J3accel( r )
            mu = 398600 ;
            J3 = 2.53241051856772e-6 ;
            Re = 6378 ;
            aP = zeros(3,1) ;
            rmag = norm(r) ;
            aP(1) = -5*J3*mu*Re^3*r(1)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
            aP(2) = -5*J3*mu*Re^3*r(2)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
            aP(3) = -5*J3*mu*Re^3/(2*rmag^7)*(6*r(3)^2-7*r(3)^4/rmag^2 - (3/5)*rmag^2) ;
        end
end