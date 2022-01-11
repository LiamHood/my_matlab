function [ rho , az , el , drho , daz , del , VIS ] = RAZEL( r , v , UTC , lat , long , alt ) 
% find azimuth and elevation
    Re = 6378 ;
    d2r = pi/180 ;
    JD = juliandate( UTC ) ;
    dcm = dcmeci2ecef( 'IAU-2000/2006' , UTC ) ;
    eci2ecef = dcm ;
    recef = dcm*r ;
    vecef = dcm*v ;
    rsite_ecef = ( lla2ecef( [ lat , long , alt ] )*1e-3 )' ;
    rhoecef = recef - rsite_ecef ;
    drhoecef = vecef ;
    ecef2sez = [ sin(lat*d2r)*cos(long*d2r) , sin(lat*d2r)*sin(long*d2r) , -cos(lat*d2r) ; ...
                -sin(long*d2r) , cos(long*d2r) , 0 ; ...
                cos(lat*d2r)*cos(long*d2r) , cos(lat*d2r)*sin(long*d2r) , sin(lat*d2r) ] ;
    rhosez = ecef2sez*rhoecef ;
%     rhosez = angle2dcm( 0 , (90-lat)*d2r , 0 )*angle2dcm( 0 , 0 , long )*rhoecef ;
%     rhosez = angle2dcm( 0 , lat , 0 )*angle2dcm( 0 , 0 , long )*rhoecef ;
    drhosez = ecef2sez*drhoecef ;
%                 drhosez = angle2dcm( 0 , (90-lat)*d2r , 0 )*angle2dcm( 0 , 0 , long )*drhoecef ;
%                 drhosez = angle2dcm( 0 , lat , 0 )*angle2dcm( 0 , 0 , long )*drhoecef ;
                rho = norm( rhosez ) ;
                el = asin( rhosez(3)/rho ) ;
                if el ~= pi/2
                    sin_az = rhosez(2) / sqrt( rhosez(1)^2 + rhosez(2)^2 ) ;
                    cos_az = -rhosez(1) / sqrt( rhosez(1)^2 + rhosez(2)^2 ) ;
                    if sin_az >= 0 && cos_az >= 0
                        az = asin( sin_az ) ;
                    elseif sin_az >= 0 && cos_az <= 0
                        az = pi - asin( sin_az ) ;
                    elseif sin_az <= 0 && cos_az <= 0
                        az = pi - asin( sin_az ) ;
                    elseif sin_az <= 0 && cos_az >= 0
                        az = asin( sin_az ) + 2*pi ;
                    end
                else
                    sin_az = drhosez(2) / sqrt( drhosez(1)^2 + drhosez(2)^2 ) ;
                    cos_az = drhosez(1) / sqrt( drhosez(1)^2 + drhosez(2)^2 ) ;
                    if sin_az >= 0 && cos_az >= 0
                        az = asin( sin_az ) ;
                    elseif sin_az >= 0 && cos_az <= 0
                        az = pi - asin( sin_az ) ;
                    elseif sin_az <= 0 && cos_az <= 0
                        az = pi - asin( sin_az ) ;
                    elseif sin_az <= 0 && cos_az >= 0
                        az = asin( sin_az ) + 2*pi ;
                    end
                end
                drho = dot( rhosez , drhosez )/rho ;
                daz = ( drhosez(1)*rhosez(2) - drhosez(2)*rhosez(1) )/( rhosez(1)^2 + rhosez(2)^2 ) ;
                del = ( rhosez(3) - drho*sin(el) )/sqrt( rhosez(1)^2 + rhosez(2)^2 ) ;

    if rhosez(3) >= 0
        rs = SunVector( JD ) ;
        if dot( rs , rsite_ecef ) > 0 
            VIS = "Radar Sun" ;
        else 
            ang = asin( norm( cross( rs , recef ) )/( norm( rs )*norm( recef ) ) ) ;
            dist = norm( recef )*cos( ang - pi/2 ) ;
            if dist > Re
                VIS = "Visible" ;
            else
                VIS = "Radar Night" ;
%                 VIS = "Visible" ;
            end
        end
    else
        VIS = "Obscured" ;
    end
end