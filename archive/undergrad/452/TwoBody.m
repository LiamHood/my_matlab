function [ t , tadj , rv , COES ] = TwoBody( tspan , r0v0 , mu , tol , sf )
fprintf( "Starting Two Body \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol ) ;
    [ t , rv ] = ode45( @TwoBodyForceFun , tspan , r0v0 , opts , mu ) ;
    ii = 2 ;
    tadj(1:2) = t(1:2) ;
    COES(1,:) = state2COE( rv(1,:) , mu ) ;
    sf(1) = floor(length(t)/1e4) ;
    sf(2) = floor(length(t)/1e4) ;
    sf(3) = floor(length(t)/1e4) ;
    while tadj(ii-1) < t( floor( end/2 - sf(1) ) ) 
        index1 = (ii-1)*sf(1) ;
        COES(ii,:) = state2COE( rv(index1,:) , mu ) ;
        tadj(ii) = t( index1 ) ;
        ii = ii + 1 ;
    end
    jj = 0 ;
    while tadj(ii + jj-1) < t( floor( end*3/4 - sf(2) ) ) 
        index2 = index1 + (jj)*sf(2);
        COES(ii+jj,:) = state2COE( rv( index2,:) , mu ) ;
        tadj(ii+jj) = t( index2 ) ;
        jj = jj + 1 ;
    end
    jj = jj + ii ;
    kk = 0 ;
    while tadj( jj + kk-1 ) < t( floor( end - sf(3) ) ) 
        index3 = index2 + (kk)*sf(3);
        COES(jj+kk,:) = state2COE( rv( index3,:) , mu ) ;
        tadj(jj+kk) = t( index3 ) ;
        kk = kk + 1 ;
    end
    
    
    
    function drdv = TwoBodyForceFun( t , rv , mu )
        r = rv(1:3) ;
        v = rv(4:6) ;

        dr = v ;
        dv = ( -mu / norm(r)^3 )*r ;
        drdv = [ dr ; dv ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , rv , mu )
        r = norm( rv(1:3) ) ;
        value = r - ( 1737.4 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end

end