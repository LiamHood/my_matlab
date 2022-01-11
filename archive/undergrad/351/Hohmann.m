function [ h , ecc , va , vp ] = Hohmann( mu , ra , rp )

    ecc = ( ra - rp ) / (  ra + rp ) ;
    h = sqrt( rp*mu*(1+ecc) ) ;
    
    va = (mu/h)*(1-ecc) ;
    vp = (mu/h)*(1+ecc) ;
    
end