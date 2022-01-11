function [ v1 , v2 ] = LambertsBook( r1 , r2 , dt , mu , tol , pro )
    
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
  
    
    A = sin( dtheta )*sqrt( r1mag*r2mag/(1-cos(dtheta)) ) ;
       
    z = -100 ;
    while F(z,dt) < 0
        z = z + .1 ;
    end
    nmax = 5000 ;
    ratio = 1 ;
    n = 0 ;
    while ( abs(ratio) > tol ) & ( n <= nmax )
        n = n + 1 ;
        ratio = F(z,dt)/dFdz(z) ;
        z = z - ratio ;
    end
        
        f = 1 - y(z)/r1mag ;
        g = A*sqrt(y(z)/mu) ;
        gdot = 1 - y(z)/r2mag ;
        
        v1 = ( 1/g )*( r2 - f*r1 ) ;
        v2 = ( 1/g )*( gdot*r2 - r1 ) ;

    
    function dum = y(z)
        [ S , C ] = Stumpf( z ) ;
        dum = r1mag + r2mag + A*(z*S-1)/sqrt(C) ;
    end
    
    function dum = F(z,dt)
        [ S , C ] = Stumpf( z ) ;
        dum = (y(z)/C)^1.5*S + A*sqrt(y(z)) - sqrt(mu)*dt ;
    end

    function dum = dFdz(z)
        [ S , C ] = Stumpf( z ) ;
        if z == 0 
            dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0))) ;
        else
            dum = (y(z)/C)^1.5*(1/2/z*(C-3*S/2/C) + 3*S^2/4/C) + A/8*(3*S/C*sqrt(y(z) + A*sqrt(C/y(z)))) ;
        end
    end
   function [ S , C ] = Stumpf( z ) 
% Stumpff Functions
sl = 11 ;
sc = zeros(1,sl) ;
selement = zeros(1,sl) ;
cc = zeros(1,sl) ;
celement = zeros(1,sl) ;
S = 0 ;
C = 0 ;

        for kk = 1:sl
            sc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+3) ;
        end
        for kk = 1:sl
            cc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+2) ;
        end
        
        for jj = 1:sl
            selement(jj) = sc(jj)*z^(jj-1) ;
        end
        S = sum(selement) ;
        
        for jj = 1:sl
            celement(jj) = cc(jj)*z^(jj-1) ;
        end
        C = sum(celement) ;
end
end