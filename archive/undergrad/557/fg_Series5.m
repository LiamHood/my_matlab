function [ f , g ] = fg_Series5( r , v , mu , t )
    rmag = norm(r) ;
    vmag = norm(v) ;
        u = mu/rmag^3 ;
    a = -u*r ;
    amag = norm(a) ;
        usd = ( 3*mu*vmag )/rmag^4 ;
    da = -usd*r - u*v ;
    damag = norm(da) ;
        udd = 3*mu*( rmag*amag - 4*vmag^2 )/( rmag^5 ) ;
        utd = 3*mu*( rmag^2*damag + 20*vmag^3 - 12*rmag*vmag*amag )/rmag^6 ;
    
    f1 = 1 ;
    f2 = ( u/2 )*t^2 ;
    f3 = ( usd/6 )*t^3 ;
    f4 = ( ( udd - u^2 )/ 20 )*t^4 ;
    f5 = ( ( utd - u^2 )/ 120 )*t^5 ;
    f = f1 - f2 - f3 - f4 - f5 ;
    
    g1 = t ;
    g2 = ( u/6 )*t^3 ;
    g3 = ( usd/12 )*t^4 ;
    g4 = ( ( 3*udd - u^2 )/ 120 )*t^5 ;
    g5 = ( ( 2*utd - 3*u*usd )/ 360 )*t^6 ;
    g = g1 - g2 - g3 - g4 - g5 ;

end