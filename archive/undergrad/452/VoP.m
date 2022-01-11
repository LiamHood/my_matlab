function [ t , COES ] = VoP( tspan , r0v0 , mu , tol , forces , A , m , epoch , JDo ) 
    COESo = state2COE( r0v0 , mu ) ;
    fprintf( "Starting VoP \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    [ t , COES ] = ode45( @VoPpropagation , tspan , COESo , opts , mu , forces , A , m , epoch , JDo) ;
    COES(:,7) = ( COES(:,1).^2 ./ mu ).*( 1 ./ ( 1 - COES(:,3).^2 ) ) ;
    COES(:,8) = COES(:,7).*(1-COES(:,3)) ;
    COES(:,9) = COES(:,7).*(1+COES(:,3)) ;
    
    function dCOES = VoPpropagation( t , COES , mu , forces , A , m , epoch , JDo )
        
        
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
        
        apeci = PerturbedAccelerations( forces , r , v , A , m , JDo , t )  ;
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

    function [ value , isterminal , direction ] = Reentry( t , COES , mu , forces , A , m , epoch , JDo )
        a = ( COES(1).^2 ./ mu ).*( 1 ./ ( 1 - COES(3).^2 ) ) ;
        rp = a.*(1-COES(3)) ;
        value = rp - ( 6378 + 100 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end
end