function [t, rv, COES] = VoP(tspan, r0v0, mu, tol, forces, A, m, Re) 
    COESo = state2coes( r0v0 , mu ) ;
    fprintf( "Starting VoP \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    [ t , COES ] = ode45( @VoPpropagation , tspan , COESo , opts , mu , forces , A , m , Re) ;
    COES(:,7) = ( COES(:,1).^2 ./ mu ).*( 1 ./ ( 1 - COES(:,3).^2 ) ) ;
    COES(:,8) = COES(:,7).*(1-COES(:,3)) ;
    COES(:,9) = COES(:,7).*(1+COES(:,3)) ;
    
    for ii = 1:length(t)
        rv(:,ii) = coes2state(COES(ii,:), mu);
    end
    function dCOES = VoPpropagation( t , COES , mu , forces , A , m , Re )
        t/(3600*24*365)
        h = COES(1) ;
        inc = COES(2) ;
        ecc = COES(3) ;
        raan = COES(4) ;
        omega = COES(5) ;
        theta = COES(6) ;
        a = COES(7) ;
        rp = COES(8) ;
        ra = COES(9) ;
        
        [ r , v ] = coes2state( COES , mu ) ;
        Nhat = cross( r , v )/norm( cross( r , v ) ) ;
        Rhat = r/norm(r) ;
        That = cross( Nhat , Rhat )/norm( cross( Nhat , Rhat ) ) ;
        Crtn = [ Rhat , That , Nhat ] ;
            if contains( forces , "drag" )
                apd = DragAcceleration( r , v , A , m ) ;
            else
                apd = zeros(3,1) ;
            end
            if contains( forces , "gravity" )
                if contains( forces , "J2" )
                    aph = J2accel( r ) ; 
                end
                if contains( forces , "J3" )
                    aph = aph + J3accel( r ) ; 
                end
            else 
                aph = zeros(3,1) ; 
            end
            apeci = apd + aph;
            ap = Crtn'*apeci ;
            R = ap(1) ;
            T = ap(2) ;
            N = ap(3) ;
            rmag = norm(r) ;
            u = omega + theta ;
            
            dh = rmag*T ;
            decc = ( h/mu ).*R.*sin(theta) + (1/(mu*h))*((h^2+mu*rmag)*cos(theta) + mu*ecc*rmag )*T ;
            dthetap = (1/(ecc*h))*( (h^2/mu)*R*cos(theta) - ((h^2/mu)+rmag)*T*sin(theta) ) ;
            dtheta = h/rmag^2 + dthetap ;
            dinc = (rmag/h)*N*cos(u) ;
            draan = ( (rmag*sin(u))/(h*sin(inc)) )*N ;
            domega = -rmag*sin(u)/(h*tan(inc))*N - dthetap ;
            
            da = 0 ;
            drp = 0 ;
            dra = 0 ;
            dCOES = [ dh ; dinc ; decc ; draan ; domega ; dtheta ; da; drp ; dra ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , COES , mu , forces , A , m , Re )
        a = ( COES(1).^2 ./ mu ).*( 1 ./ ( 1 - COES(3).^2 ) ) ;
        rp = a.*(1-COES(3)) ;
        value = rp - ( Re + 120 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end
end