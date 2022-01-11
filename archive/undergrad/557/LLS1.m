function [ yc , rbar , RMS , P ] = LLS1( xoi , yoi )
        n = length( xoi ) ;
        AtA = zeros( 2 ) ;
        AtA(1,1) = n ;
        AtA(1,2) = sum( xoi ) ;
        AtA(2,1) = sum( xoi ) ;
        AtA(2,2) = sum( xoi.^2 ) ;
        Atb = zeros( 2 , 1 ) ;
        Atb(1,1) = sum( yoi ) ;
        Atb(2,1) = sum( yoi.*xoi ) ;
        P = inv( AtA ) ;
        state = P*Atb ;
        alpha = state(1) ;
        beta = state(2) ;
        yc = alpha + beta*xoi ;
        rbar = yoi - yc ;
        RMS = sqrt( sum( rbar.^2 )/n ) ;
end