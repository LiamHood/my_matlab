function [ vp , va , vinner , vouter , dvi , dvo ] = HohmannC2C( mu , router , rinner )

    ecc = ( router - rinner ) / (  router + rinner ) ; 
    h = sqrt( rinner*mu*(1+ecc) ) ; 
    
    va = (mu/h)*(1-ecc) ;
    vp = (mu/h)*(1+ecc) ;
    
    vinner = sqrt( mu / rinner ) ;
    vouter = sqrt( mu / router ) ;
    
    dvi = vp - vinner ;
    dvo = vouter - va ;
    
end