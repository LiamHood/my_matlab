function [ r , v ] = Pcoes2state( p_coes , mu )
        d2r = pi/180 ;
        h = sqrt( p_coes(1)*mu*(1-p_coes(2)^2) ) ;
        omega = p_coes(6) - p_coes(4) ;
        M = p_coes(6)*d2r - p_coes(5)*d2r ;
        T = (2*pi)/sqrt(mu)*p_coes(1) ;
        t = M*T/(2*pi) ;
        [ theta ] = time2theta( t , T , p_coes(2) ) ;
        [ r , v ] = coes2state(  h , p_coes(2) , theta , p_coes(4) , omega , p_coes(3) , mu ) ;
                
    function [ theta ] = time2theta( t , T , ecc )
    % Find true anomaly at a time

    n = 2*pi/T ; % mean motion
    Me = n*t ; 

    % Guess of E
    if Me < pi
        E0 = Me + ecc/2 ;
    else
        E0 = Me - ecc/2 ;
    end

    % Use Newtons to find E
        tol = 10^-8 ; % Tolerance
        lim = 1000 ; % Maximum iteration
        f = @(E) E - ecc*sin(E) - Me ; % Function handle for E
        fprime = @(E) 1 - ecc*cos(E) ; % function handle for derivative of E
        [ E ] = newton( E0 , f , fprime , tol , lim ) ; % Apply Newtons

    theta = 2*atan(tan(E/2)*sqrt((1+ecc)/(1-ecc))) ; % find true anomaly
    % correction to make it positive
        if theta < 0
            theta = theta + 2*pi ;
        end
    theta = theta*(180/pi) ;
    end
    function [ r , v ] = coes2state(  h , ecc , theta , RAAN , omega , inc , mu )
        r_peri = (h^2/mu) * ( 1/( 1 + ecc*cosd(theta) ) ) * [ cosd( theta ) ; sind( theta ) ; 0 ] ;
        v_peri = (mu/h) * [ -sind( theta ) ; ecc+cosd(theta) ; 0 ] ;

        d2r = pi/180 ;
        RAAN = d2r*RAAN ;
        omega = d2r*omega ;
        inc = d2r*inc ;
        Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
        Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
        Q(1,3) = sin(RAAN)*sin(inc) ;
        Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
        Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
        Q(2,3) = -cos(RAAN)*sin(inc) ;
        Q(3,1) = sin(inc)*sin(omega) ;
        Q(3,2) = sin(inc)*cos(omega) ;
        Q(3,3) = cos(inc) ;

        r = Q*r_peri ;
        v = Q*v_peri ;
    end
end