function [ f , g ] = fg_Series2( r , mu , t )
    rmag = norm(r) ;
        
    u = mu/rmag^3 ;
    
    f1 = 1 ;
    f2 = ( u/2 )*t^2 ;
    f = f1 - f2 ;
    
    g1 = t ;
    g2 = ( u/6 )*t^3 ;
    g = g1 - g2 ;

end