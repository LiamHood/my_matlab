function [ rmag , ra , dec , drmag , dra , ddec ] = Geocentric_RA_DEC( r , v ) 
% Finds geocentric right ascension and declination 
    rmag = norm( r ) ;
    sin_dec = ( r(3)/rmag ) ;
    cos_dec = ( sqrt( r(1)^2 + r(2)^2 )/rmag ) ;
    if sin_dec >= 0 && cos_dec >= 0
        dec = asin( sin_dec ) ;
    elseif sin_dec >= 0 && cos_dec <= 0
        dec = pi - asin( sin_dec ) ;
    elseif sin_dec <= 0 && cos_dec <= 0
        dec = pi - asin( sin_dec ) ;
    elseif sin_dec <= 0 && cos_dec >= 0
        dec = asin( sin_dec ) + 2*pi ;
    end
    if sqrt( r(1)^2 + r(2)^2 ) ~= 0
        sin_ra = r(2) / sqrt( r(1)^2 + r(2)^2 ) ;
        cos_ra = r(1) / sqrt( r(1)^2 + r(2)^2 ) ;
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
        sin_ra = v(2) / sqrt( v(1)^2 + v(2)^2 ) ;
        cos_ra = v(1) / sqrt( v(1)^2 + v(2)^2 ) ;
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
    drmag = dot( r , v )/rmag ;
    dra = ( v(1)*r(2) - v(2)*r(1) )/( - r(2)^2 - r(1)^2 ) ;
    ddec = ( v(3) - drmag*( r(3) / rmag ) )/sqrt( r(1)^2 + r(2)^2 ) ;
end