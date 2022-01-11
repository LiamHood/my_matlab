function [ r , v ] = NewState( r0 , v0 , dt , mu )
% Find new position and velocity at some time in orbit by lagrange
% variables

    % Initial values
    r0mag = norm( r0 ) ; % radius in km
    v0mag = norm( v0 ) ; % speed in km/s
    vr0 = dot( r0 , v0 ) / r0mag ; % radial speed (km/s)
    alpha = (2/r0mag) - ( v0mag^2 / mu ) ;
    chi = sqrt( mu ) * abs( alpha ) * dt ; % Initial universal anomaly guess
    z = alpha*chi(1)^2 ;
    sl = 10 ; % series length
    
    % Universal anomaly equations
    fun = @(chi,C,S) ((r0mag*vr0)/sqrt(mu))*chi^2*C+(1-alpha*r0mag)*chi^3*S+r0mag*chi-sqrt(mu)*dt ;
    fprime = @(chi,C,S) ((r0mag*vr0)/sqrt(mu))*chi*(1-alpha*chi^2*S)+(1-alpha*r0mag)*chi^2*C+r0mag ;
    
    % Stumpff Functions
        for kk = 1:sl
            sc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+3) ;
        end
        

        for kk = 1:sl
            cc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+2) ;
        end
        

    % Newtons for UV
    ii = 1 ;
    ratio = 1 ;
    tol = 10^-8 ;
    lim = 1000 ;
        while abs(ratio(ii)) >= tol
            
            % Stumpff calc
                for jj = 1:sl
                    S = sc(jj)*z^(jj-1) ;
                end
                for jj = 1:sl
                   C = cc(jj)*z^(jj-1) ;
                end
                
            ratio(ii+1) = fun(chi,C,S)/fprime(chi,C,S) ;
            chi = chi - ratio(ii+1) ;
            ii = ii + 1 ;
            z = alpha*chi^2 ; % New z
                if ii > lim
                    error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
                end
        end
    % Final Stumpff
    for jj = 1:sl
        S = sc(jj)*z^(jj-1) ;
    end
    for jj = 1:sl
        C = cc(jj)*z^(jj-1) ;
    end
    
    
    % Lagrange variables
    
    % new r
    f = 1 - (chi^2/r0mag)*C ;
    g = dt - (1/sqrt(mu))*chi^3*S ;
    r = f*r0 + g*v0 ; % new position
    rmag = norm( r ) ; % radius
    
    % new v
    fdot = (sqrt(mu)/(r0mag*rmag))*(alpha*chi^3*S-chi) ;
    gdot = 1 - (chi^2/rmag)*C ; 
    v = fdot*r0 + gdot*v0 ; % new velocity
    

end