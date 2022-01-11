function [ v2 ] = Gibbs( r1 , r2 , r3 , JD1 , JD2 , JD3 , mu ) 
        
        z12 = cross( r1 , r2 ) ;
        z23 = cross( r2 , r3 ) ;
        z31 = cross( r3 , r1 ) ;
        r1mag = norm( r1 ) ;
        r2mag = norm( r2 ) ;
        r3mag = norm( r3 ) ;

        alpha_cop = asind( dot( z23 , r1 )/( norm( z23 )*norm( r1 ) ) ) ;
            % check seperation
            alpha12 = acos( dot( r1 , r2 )/( norm( r1 )*norm( r2 ) ) ) ;
            alpha23 = acos( dot( r2 , r3 )/( norm( r2 )*norm( r3 ) ) ) ;
            ang_sep = abs( alpha12 - alpha23 )*( 180/pi ) ;
        if ang_sep > 5
            N = r1mag*z23 + r2mag*z31 + r3mag*z12 ;
            D = z12 + z23 + z31 ;
            S = ( r2mag - r3mag )*r1 + ( r3mag - r1mag )*r2 + ( r1mag - r2mag )*r3 ;
            B = cross( D , r2 ) ;
            Lg = sqrt( mu/( norm( N )*norm( D ) ) ) ;
            v2 = ( Lg / r2mag )*B + Lg*S ;
        else
            dt31 = ( JD3 - JD1 )*86400 ;
            dt32 = ( JD3 - JD2 )*86400 ;
            dt21 = ( JD2 - JD1 )*86400 ;
            alpha_cop = 90 - acosd( dot( z23 , r1 )/( norm( z23 )*norm( r1 ) ) ) ;
            v2 = -dt32*( ( 1/(dt21*dt31) ) + ( mu/(12*r1mag^3) ) )*r1 + ...
                (dt32-dt21)*( ( 1/(dt21*dt32) ) + ( mu/(12*r2mag^3) ) )*r2 + ...
                dt21*( ( 1/(dt32*dt31) ) + ( mu/(12*r3mag^3) ) )*r3 ;
        end
    end