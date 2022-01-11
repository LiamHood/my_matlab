function [ r2 , v2 , r1 , r3 ] = GaussExtended( Lhat1 , Lhat2 , Lhat3 , R1 , R2 , R3 , tau1 , tau3 ) 
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
    %-----------------------Extension------------------------------
    var = 1 ;
    tol = 1e-4 ;
    iteration = 0 ;
    lim = 1e4 ;
    rho2(2) = rho2 ;
    while var >= tol 
        rho2(1) = rho2(2) ;
        
        [ f1 , g1 ] = LagrangianUV( r2 , v2 , tau1 , mu ) ;
        [ f3 , g3 ] = LagrangianUV( r2 , v2 , tau3 , mu ) ;

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