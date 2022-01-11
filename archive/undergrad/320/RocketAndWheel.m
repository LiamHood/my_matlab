function [ dstate ] = RocketAndWheel( t , state , dwrel , Td , Iw , Ir )

    thetar = state(1:3) ;
    wr = state(4:6) ;
    wrel = state(7:9) ;
    dstate = zeros( 9 , 1 ) ;

    dwr = -inv( Ir+Iw ) * ( Iw*dwrel + crossmatrix(wr+wrel)*Iw*(wr+wrel) + crossmatrix(wr)*Ir*wr - Td ) ;
    dstate(1:3) = wr ;
    dstate(4:6) = dwr ;
    dstate(7:9) = dwrel ;
end