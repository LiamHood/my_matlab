function [ r2 , v2 , r1 , r3 ] = GaussIOD( Lhat1 , Lhat2 , Lhat3 , R1 , R2 , R3 , tau1 , tau3 ) 
    % finds state of middle observation using gauss angle only method. Lhat
    % is the slant range unit vector, tau1,3 is time1,3 difference from
    % time2 in seconds, R is the eci site vector in km
    mu = 398600 ;
    a1 = tau3/( tau3 - tau1 ) ;
    a1u = tau3*( ( tau3 - tau1 )^2 - tau3^2 )/( 6*( tau3 - tau1 ) ) ;
    a3 = -tau1/( tau3 - tau1 ) ;
    a3u = - tau1*( ( tau3 - tau1 )^2 - tau1^2 )/( 6*( tau3 - tau1 ) ) ;

    M = ( [ Lhat1 , Lhat2 , Lhat3 ] )\[ R1 , R2 , R3 ] ;
    d1 = M(2,1)*a1 - M(2,2) + M(2,3)*a3 ;
    d2 = M(2,1)*a1u + M(2,3)*a3u ;
    C = dot( Lhat2 , R2 ) ;
    
    % Find magnitude of r2 by finding real positive root
    syms x ;
    eqn = x^8 - ( d1^2 + 2*C*d1 + norm(R2)^2 )*x^6 - 2*mu*( C*d2 + d1*d2 )*x^3 - ( mu^2 * d2^2 ) == 0 ;
    y = double( subs( vpa( solve( eqn , x , 'Real' , true ) ) ) ) ;
    for ii = 1:length(y)
        if y(ii) > 0
            r2mag = y(ii) ;
        end
    end
    
    u = mu/r2mag^3 ;
    c1 = a1 + a1u*u ;
    c2 = -1 ;
    c3 = a3 + a3u*u ;
    
    crho = M*[ -c1 ; -c2 ; -c3 ] ;
    rho1 = crho(1)/c1 ;
    rho2 = crho(2)/c2 ;
    rho3 = crho(3)/c3 ;
    
    r1 = rho1*Lhat1 + R1 ;
    r2 = rho2*Lhat2 + R2 ;
    r3 = rho3*Lhat3 + R3 ;
    r2ne = r2 ;

    [ f1 , g1 ] = fg_Series2( r2 , mu , tau1 ) ;
    [ f3 , g3 ] = fg_Series2( r2 , mu , tau3 ) ;
    
    v2 = (1/( f1*g3 - f3*g1 ))*( -f3*r1 + f1*r3 ) ;
    v2ne = v2 ;

    %-------------------Gibbs Extension------------------------------------
    var = 1 ;
    tol = 1e-8 ;
    iteration = 0 ;
    lim = 1e4 ;
    rho2(2) = rho2 ;
    while var >= tol 
        rho2(1) = rho2(2) ;
        
        z12 = cross( r1 , r2 ) ;
        z23 = cross( r1 , r2 ) ;
        z31 = cross( r1 , r2 ) ;
        r1mag = norm( r1 ) ;
        r2mag = norm( r2 ) ;
        r3mag = norm( r3 ) ;

        alpha_cop = asin( dot( z23 , r1 )/( norm( z23 )*norm( r1 ) ) ) ;
            % check seperation
            alpha12 = acos( dot( r1 , r2 )/( norm( r1 )*norm( r2 ) ) ) ;
            alpha23 = acos( dot( r2 , r3 )/( norm( r2 )*norm( r3 ) ) ) ;
            ang_sep = abs( alpha12 - alpha23 )*( 180/pi ) ;
        if ang_sep > 1
            N = r1mag*z23 + r2mag*z31 + r3mag*z12 ;
            D = z12 + z23 + z31 ;
            S = ( r2mag - r3mag )*r1 + ( r3mag - r1mag )*r2 + ( r1mag - r2mag )*r3 ;
            B = cross( D , r2 ) ;
            Lg = sqrt( mu/( norm( N )*norm( D ) ) ) ;
            v2 = Lg/r2mag*B + Lg*S ;
        else
            error( 'Add Herrick Gibbs' )
        end

        p = ( norm( cross( r2 , v2 ) )^2 ) / mu ;
        f1 = 1 - ( r1mag / p )*( 1 - cos( alpha12 ) ) ;
        f3 = 1 - ( r3mag / p )*( 1 - cos( alpha23 ) ) ;
        g1 = r1mag*r2mag*sin( alpha12 )/sqrt( mu*p ) ;
        g3 = r3mag*r2mag*sin( alpha23 )/sqrt( mu*p ) ;
        c1 = g3/( f1*g3 - f3*g1 ) ;
        c2 = -1 ;
        c3 = -g1/( f1*g3 - f3*g1 ) ;
        crho = M*[ -c1 ; -c2 ; -c3 ] ;
        rho1 = crho(1)/c1 ;
        rho2(2) = crho(2)/c2 ;
        rho3 = crho(3)/c3 ;

        r1 = rho1*Lhat1 + R1 ;
        r2 = rho2(2)*Lhat2 + R2 ;
        r3 = rho3*Lhat3 + R3 ;
        v2 = 1/( f1*g3 - f3*g1 )*( -f3*r1 + f1*r3 ) ;
        
        var = abs( rho2(2) - rho2(1) ) ;
        iteration = iteration + 1 ;
                    if iteration > lim
                        error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
                    end
    end
 
end