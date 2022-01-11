function [ rho , ra , dec , drho , dra , ddec ] = Topocentric_RA_DEC( r , v , rsite , vsite )
% Finds topocentric right ascension and declination 
    rhov = r - rsite ;
    rho = norm( rhov ) ;
    dec = asin( rhov( 3 )/rho ) ;
    drhov = v - vsite ;
    if sqrt( rhov(1)^2 + rhov(2)^2 ) ~= 0
        sin_ra = rhov(2) / sqrt( rhov(1)^2 + rhov(2)^2 ) ;
        cos_ra = rhov(1) / sqrt( rhov(1)^2 + rhov(2)^2 ) ;
        if sin_ra >= 0 && cos_ra >= 0
            ra = asin( sin_ra ) ;
        elseif sin_ra >= 0 && cos_ra <= 0
            ra = pi - asin( sin_ra ) ;
        elseif sin_ra <= 0 && cos_ra <= 0
            ra = pi - asin( sin_ra ) ;
        elseif sin_ra <= 0 && cos_ra >= 0
            ra = asin( sin_ra ) + 2*pi ;
        end
    else
        sin_ra = drhov(2) / sqrt( drhov(1)^2 + drhov(2)^2 ) ;
        cos_ra = drhov(1) / sqrt( drhov(1)^2 + drhov(2)^2 ) ;
        if sin_ra >= 0 && cos_ra >= 0
            ra = asin( sin_ra ) ;
        elseif sin_ra >= 0 && cos_ra <= 0
            ra = pi - asin( sin_ra ) ;
        elseif sin_ra <= 0 && cos_ra <= 0
            ra = pi - asin( sin_ra ) ;
        elseif sin_ra <= 0 && cos_ra >= 0
            ra = asin( sin_ra ) + 2*pi ;
        end
    end
    drho = dot(  rhov , drhov )/rho ;
    dra = ( drhov(1)*rhov(2) - drhov(2)*rhov(1) )/( - rhov(2)^2 - rhov(1)^2 ) ;
    ddec = ( drhov(3) - drho*sin( dec ) )/sqrt( rhov(1)^2 + rhov(2)^2 ) ;

end