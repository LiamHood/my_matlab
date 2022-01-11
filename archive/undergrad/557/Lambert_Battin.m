function [ v1 , v2 ] = Lambert_Battin( r1 , r2 , dt , tm , mu )
% uses universal variable to find velocity at two postions 
% r1,2 are position vectors, columns, in km
% dt is time of flight in seconds
% tm is short or long way, 1 for short or -1 for long
% v1,v2 are column vectors in km/s
    r1mag = norm( r1 ) ;
    r2mag = norm( r2 ) ;
    cos_ta = dot( r1 , r2 )/( r1mag*r2mag ) ;
    sin_ta = tm*sqrt( 1 - cos_ta^2 ) ;
    if sin_ta > 0
        ta = acos( cos_ta ) ;
    else
        ta = 2*pi + asin( sin_ta ) ;
    end
    c = sqrt( r1mag^2 + r2mag^2 - 2*r1mag*r2mag*cos_ta ) ;
    s = ( r1mag + r2mag + c )/2 ;
    eps = ( r2mag - r1mag )/r1mag ;
    tan2_2w = ( eps^2 / 4 )/( sqrt( r2mag/r1mag ) + ( r2mag/r1mag )*( 2 + ( sqrt( r2mag/r1mag ) ) ) ) ;
    rop = sqrt( r1mag*r2mag )*( cos( ta/4 )^2 + tan2_2w ) ;
    if ta < 180 
        L = ( sin( ta/4 )^2 + tan2_2w )/( sin( ta/4 )^2 + tan2_2w + cos( ta/2 ) ) ;
    else 
        L = ( cos( ta/4 )^2 + tan2_2w - cos( ta/2 ) )/( cos( ta/4 )^2 + tan2_2w ) ;
    end
    m = ( mu*dt^2 )/( 8*rop^3 ) ;
    x = L ;
    xv = [ L , L ] ;
    var = 1 ;
    tol = 1e-6 ;
    iteration = 0 ;
    lim = 1e3 ;
    while var > tol
        xv(1) = xv(2) ;
        x = xv(2) ;
        eta = x/( ( sqrt( 1 + x ) + 1 )^2 ) ;
        ksi = 8*( sqrt( 1 + x ) + 1 )/( 3 + ( 1/( 5 + eta + ( (9/7)*eta/...
            (1+((16/63)*eta/(1+((25/99)*eta/(1+((36/143)*eta/(1+((49/195)*eta/...
            (1+((64/255)*eta/1.25)))))))))))))) ;
        h1 = ( ( L + x )^2 * ( 1 + 3*x + ksi ) )/( ( 1 + 2*x + L )*( 4*x + ksi*( 3 + x ) ) ) ;
        h2 = ( m*( x - L + ksi ) )/( ( 1 + 2*x + L )*( 4*x + ksi*( 3 + x ) ) ) ;
        B = ( 27*h2 )/( 4*( 1 + h1 )^3 ) ;
        U = B/( 2*( sqrt( 1 + B ) + 1 ) ) ;
        K = ( 1/3 )/(1+((4/27)*U)/(1+((8/27)*U)/(1+(0.2334*U)/(1+(0.2642*U)/(1+(0.2408*U)/(1+(0.2584*U)/(1+(0.2436*U)/1.25))))))) ;
        y = ( ( 1 + h1 )/3 )*( 2 + sqrt( 1 + B )/( 1 + 2*U*K^2 ) ) ;
        x = sqrt( ( ( 1 - L )/2 )^2 + ( m / y^2 ) ) - ( 1 + L )/2 ;
%         syms cubic ;
%         eqn = cubic^3 - cubic^2 - h1*cubic^2 - h2 == 0 ;
%         yg = double( subs( vpa( solve( eqn , cubic , 'MaxDegree' , 3 , 'Real' , true ) ) ) ) ;
%         for ii = 1:length(yg)
%             if yg(ii) > 0 
%                 y = yg(ii) ;
%             end
%         end
%         x = sqrt( ( ( 1 - L )/2 )^2 + ( m/y^2 ) ) - ( 1 + L )/2 ;
        xv(2) = x ;
        var = abs( xv(2) - xv(1) ) ;
        if iteration > lim
            error( 'dt value failed to converge' )
        end
    end
    a = ( mu*dt^2 )/( 16*rop^2*x*y^2 ) ;
    if a > 1e-12
        beta = 2*asin( sqrt( ( s - c )/( 2*a ) ) ) ;
        if ta > pi 
            beta = -beta ;
        end
        amin = s/2 ;
        tmin = sqrt( amin^3 / mu )*( pi - beta + sin( beta ) ) ;
        alpha = 2*asin( sqrt( s/(2*a) ) ) ;
        if dt > tmin
            alpha = 2*pi - alpha ;
        end
        dE = alpha - beta ;
        f = 1 - ( a/r1mag )*( 1 - cos( dE ) ) ;
        g = dt - sqrt( a^3 / mu )*( dE - sin( dE ) ) ;
        gdot = 1 - ( a/r2mag )*( 1 - cos( dE ) ) ;
    else 
        alphah = 2*asinh( sqrt( s/( -2*a ) ) ) ;
        betah = 2*asinh( sqrt( ( s - c )/( -2*a ) ) ) ;
        dH = alphah - betah ;
        f = 1 - ( a/r1mag )*( 1 - cosh( dH ) ) ;
        g = dt - sqrt( (-a)^3 / mu )*( sinh( dH ) - dH ) ;
        gdot = 1 - ( a/r2mag )*( 1 - cosh( dH ) ) ;
    end
    v1 = ( r2 - f*r1 )/g ;
    v2 = ( gdot*r2 - r1 )/g ;
end