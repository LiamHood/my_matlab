function [ r , v ] = LowThrustCanonical( dt , r0 , v0 , mu , athrust )

    tol = 1e-8 ;
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol ) ;
    [ t , rv ] = ode45( @LowThrustForceFun , [0,dt] , [ r0 ; v0 ] , opts , mu , athrust) ;
    r = rv(:,1:3)' ;
    v = rv(:,4:6)' ;
    
    function drdv = LowThrustForceFun( t , rv , mu , athrust )
        rvec = rv(1:3) ;
        vvec = rv(4:6) ;

        dr = vvec ;
        dv = ( -mu / norm(rvec)^3 )*rvec + athrust ;
        drdv = [ dr ; dv ] ;
    end


                    
end