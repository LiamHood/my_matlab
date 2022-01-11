function [ v1_long , v2_long , v1_short , v2_short ] = Lamberts( r1 , r2 , dt , mu , tol , pro )
% pro is 1 or 0 for prograde or retrograde respectively
    
    r1mag = norm( r1 ) ;
    r2mag = norm( r2 ) ;
    rcross = cross( r1 , r2 ) ;
    
    % Find delta theta
        if pro == 1 
            if rcross(3) >= 0
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        else
            if rcross(3) < 0 
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        end  
    tm = [ -1 , 1 ] ;
    for ii = 1:2
    dtheta = asin( tm(ii)*(1-cos(dtheta)^2) ) ;
    A = sin( dtheta )*sqrt( r1mag*r2mag/(1-cos(dtheta)) ) ;
        z = 0 ;
        C = 1/2 ;
        S = 1/6 ;
        zup = 4*pi^2 ;
        zlow = -4*pi^2 ;
        y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
        chi = sqrt(y/C) ;
        dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        while abs( dtloop - dt ) > tol
                if dtloop <= dt
                    zlow = z ;
                else 
                    zup = z ;
                end
             z = ( zup + zlow ) / 2 ;
            [ S , C ] = Stumpf( z ) ;
            y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
            chi = sqrt(y/C) ;
            dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        end
        f = 1 - y/r1mag ;
        g = A*sqrt(y/mu) ;
        gdot = 1 - y/r2mag ;
        
        v1(:,ii) = ( 1/g )*( r2 - f*r1 ) ;
        v2(:,ii) = ( 1/g )*( gdot*r2 - r1 ) ;
    end
    v1_long = v1(:,1) ;
    v2_long = v2(:,1) ;
    v1_short = v1(:,2) ;
    v2_short = v2(:,2) ;
end