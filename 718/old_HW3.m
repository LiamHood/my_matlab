%% Homework 3
% Aero 452
% Liam Hood
function HW3()
clear ; close all ; clc ;
mu = 398600 ;
re = 6378 ;
epoch = [ 2019 1 1 12 0 0 ] ;
pt = '\n\nProblem number %u \n \n' ;

%% Problem 1
fprintf( pt , 1 ) 
HW3_P1( mu , re , epoch )
 
%% Problem 2
fprintf( pt , 2 ) 
% HW3_P2( mu , re , epoch )

% %% Problem 3
% fprintf( pt , 3 )
% HW3_P3( mu , re , epoch )

%% Problem Solutions
 
    function HW3_P1( mu , re , epoch )
        s2d = 1/(3600*24) ;
        r2d = 180/pi ;
        d2r = pi/180 ;
        rp0 = 215 + re ;
        ra0 = 939 + re ;
        raan0 = 340*d2r ;
        inc0 = 65.2*d2r ;
        omega0 = 58*d2r ;
        theta0 = 332*d2r ;
        a0 = ( rp0 + ra0 ) / 2 ;
        ecc0 = 1 - rp0/a0 ;
        p = a0*( 1 - ecc0^2 ) ;
        h0 = sqrt( mu*p ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        COESo = [ h0 , inc0 , ecc0 , raan0 , omega0 , theta0 ] ;
        [ r0 , v0 ] = coes2state(  COESo , mu ) ;
        r0v0 = [ r0 ; v0 ] ;
        A = 1 ;
        m = 100 ;
        
        tend = 12*24*60*60 ;
        tspan = [ 0 , tend ] ;
        tol = 1e-8 ;
        dt = 50 ;
        forces = "drag" ;
        name = [ "Two Body" , "Encke" , "Cowell" , "VoP" ] ;
        sf = [ 1000 , 1000 , 100 ] ;
        sfe = [ 10 , 10 , 1 ] ;
        
        tic ;
        [ tfull2b , t2b , rv2b , COE2b] = TwoBody( tspan , r0v0 , mu , tol , sf ) ;
        t2end(1) = toc ;
        tic
        [ tfullencke , tencke , rencke , vencke , COEencke ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m , sfe ) ;
        t2end(2) = toc ;
        tic ;
        [ tfullcow , tcow , rvcow , COEcow ] = Cowell( tspan , r0v0 , mu , tol , forces , A , m , sf ) ;
        t2end(3) = toc ;
        tic ;
        [ tvop , COEvop ] = VoP( tspan , r0v0 , mu , tol , forces , A , m , epoch ) ;
        t2end(4) = toc ;
        for ii = 1:4
            fprintf( 'The time to complete %s propagation was %f seconds \n' , name(ii) , t2end(ii) )
        end
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        
        figure
        subplot( 2 , 2 , 1 )
        plot( t2b*s2d , COE2b( :,2 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        title( 'Two Body' )
        subplot( 2 , 2 , 2 )
        plot( tencke*s2d , COEencke( :,2 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        title( 'Encke' )
        subplot( 2 , 2 , 3 )
        plot( tcow*s2d , COEcow( :,2 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        title( 'Cowell' )
        subplot( 2 , 2 , 4 )
        plot( tvop*s2d , COEvop( :,2 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        title( 'Variation of Parameters' )
        
        figure
        subplot( 2 , 2 , 1 )
        plot( t2b*s2d , COE2b( :,3 ) )
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        title( 'Two Body' )
        subplot( 2 , 2 , 2 )
        plot( tencke*s2d , COEencke( :,3 ) )
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        title( 'Encke' )
        subplot( 2 , 2 , 3 )
        plot( tcow*s2d , COEcow( :,3 ) )
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        title( 'Cowell' )
        subplot( 2 , 2 , 4 )
        plot( tvop*s2d , COEvop( :,3 ) )
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        title( 'Variationof Parameters' )
        
        figure
        subplot( 2 , 2 , 1 )
        plot( t2b*s2d , COE2b( :,4 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'RAAN [degrees]'  )
        title( 'Two Body' )
        subplot( 2 , 2 , 2 )
        plot( tencke*s2d , COEencke( :,4 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'RAAN [degrees]'  )
        title( 'Encke' )
        subplot( 2 , 2 , 3 )
        plot( tcow*s2d , COEcow( :,4 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'RAAN [degrees]'  )
        title( 'Cowell' )
        subplot( 2 , 2 , 4 )
        plot( tvop*s2d , COEvop( :,4 )*r2d ) 
        xlabel( 'Time [days]' )
        ylabel( 'RAAN [degrees]' )
        title( 'Variation of Parameters' )
        
        figure
        subplot( 2 , 2 , 1 )
        plot( t2b*s2d , COE2b( :,5 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Argument of Perigee [degrees]'  )
        title( 'Two Body' )
        subplot( 2 , 2 , 2 )
        plot( tencke*s2d , COEencke( :,5 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Argument of Perigee [degrees]'  )
        title( 'Encke' )
        subplot( 2 , 2 , 3 )
        plot( tcow*s2d , COEcow( :,5 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Argument of Perigee [degrees]'  )
        title( 'Cowell' )
        subplot( 2 , 2 , 4 )
        plot( tvop*s2d , COEvop( :,5 )*r2d )
        xlabel( 'Time [days]' )
        ylabel( 'Argument of Perigee [degrees]' )
        title( 'Variation of Parameters' )
        
        figure
        subplot( 2 , 2 , 1 )
        hold on 
        plot( t2b*s2d , COE2b( :,8 )-re )
        plot( t2b*s2d , COE2b( :,9 )-re )
        ylabel( 'Altitude [km]' )
        xlabel( 'Time [days]' )
        title( 'Two Body' )
        legend( 'Perigee' , 'Apogee' )
        hold off
        subplot( 2 , 2 , 2 )
        hold on
        plot( tencke*s2d , COEencke( :,8 )-re )
        plot( tencke*s2d , COEencke( :,9 )-re )
        ylabel( 'Altitude [km]' )
        xlabel( 'Time [days]' )
        title( 'Encke' )
        legend( 'Perigee' , 'Apogee' )
        hold off
        subplot( 2 , 2 , 3 )
        hold on
        plot( tcow*s2d , COEcow( :,8 )-re )
        plot( tcow*s2d , COEcow( :,9 )-re )
        ylabel( 'Altitude [km]' )
        xlabel( 'Time [days]' )
        title( 'Cowell' )
        legend( 'Perigee' , 'Apogee' )
        hold off
        subplot( 2 , 2 , 4 )
        hold on
        plot( tvop*s2d , COEvop( :,8 )-re )    
        plot( tvop*s2d , COEvop( :,9 )-re )
        ylabel( 'Altitude [km]' )
        xlabel( 'Time [days]' )
        title( 'Variation of Parameters' )
        legend( 'Perigee' , 'Apogee' )
        hold off
        
        check = 'Each propagation technique gives similar results. This makes sense \nbecause they are supposed to be modeling the same thing \n\n' ;
        fprintf( check )
    end

   function HW3_P2( mu , re , epoch )
        s2d = 1/(3600*24) ;
        d2r = pi/180 ;
        r2d = 180/pi ;
        rp0 = 300 + re ;
        ra0 = 3092 + re ;
        raan0 = 45*d2r ;
        inc0 = 28*d2r ;
        omega0 = 30*d2r ;
        theta0 = 40*d2r ;
        a0 = ( rp0 + ra0 ) / 2 ;
        ecc0 = 1 - rp0/a0 ;
        p = a0*( 1 - ecc0^2 ) ;
        h0 = sqrt( mu*p ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        COESo = [ h0 , inc0 , ecc0 , raan0 , omega0 , theta0 ] ;
        [ r0 , v0 ] = coes2state(  COESo , mu ) ;
        r0v0 = [ r0 ; v0 ] ;
        A = 1 ;
        m = 100 ;
        
        tend = 24*60*60 ;
        tspan = [ 0 , tend ] ;
        tol = 1e-12 ;
        forces = 'gravity,J2' ;
        forces3 = 'gravity,J2,J3' ;
        sf = [ 1 , 1 , 1 ] ;
        sfe = [ 1 , 1 , 1 ] ;
        
        [ t2 , COE2 ] = VoP( tspan , r0v0 , mu , tol , forces , A , m , epoch ) ;
        [ t3 , COE3 ] = VoP( tspan , r0v0 , mu , tol , forces3 , A , m , epoch ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,2 )*r2d )
        plot( t3*s2d , COE3( :,2 )*r2d )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        legend( 'J2' , 'J2 and J3' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,3 ) )
        plot( t3*s2d , COE3( :,3 ) )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        legend( 'J2' , 'J2 and J3' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,4 )*r2d )
        plot( t3*s2d , COE3( :,4 )*r2d )
        hold off
        title( 'RAAN' )
        legend( 'J2' , 'J2 and J3' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,5 )*r2d )
        plot( t3*s2d , COE3( :,5 )*r2d )
        hold off
        title( 'Argument of Perigee [degrees]' )
        legend( 'J2' , 'J2 and J3' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,8 )-re , 'b' )
        plot( t3*s2d , COE3( :,8 )-re ,'r' )
        plot( t2*s2d , COE2( :,9 )-re , 'b' )
        plot( t3*s2d , COE3( :,9 )-re , 'r')
        hold off
        title( 'Altitude [km]' )
        legend( 'J2' , 'J2 and J3' )
        
        check = 'There is barely any change during 1 day. This makes sense because \nJ3 is orders of magnitude less than J2 \n' ;
        fprintf( check )
   end

    function HW3_P3( mu , re , epoch )
        s2d = 1/(3600*24) ;
        r2d = 180/pi ;
        d2r = pi/180 ;
        rp0 = 215 + re ;
        ra0 = 939 + re ;
        raan0 = 340*d2r ;
        inc0 = 65.2*d2r ;
        omega0 = 58*d2r ;
        theta0 = 332*d2r ;
        a0 = ( rp0 + ra0 ) / 2 ;
        ecc0 = 1 - rp0/a0 ;
        p = a0*( 1 - ecc0^2 ) ;
        h0 = sqrt( mu*p ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        COESo = [ h0 , inc0 , ecc0 , raan0 , omega0 , theta0 ] ;
        [ r0 , v0 ] = coes2state(  COESo , mu ) ;
        r0v0 = [ r0 ; v0 ] ;
        A = 1 ;
        m = 100 ;
        
        tend = 120*24*60*60 ;
        tspan = [ 0 , tend ] ;
        tol = 1e-8 ;
        forces = "drag" ;
        forces3 = "drag,NRLMSISE" ;
        name = [ "Exponential" , "NRLMSISE" ] ;
        
        tic ;
        [ t2 , COE2 ] = VoP( tspan , r0v0 , mu , tol , forces , A , m , epoch ) ;
        texp = toc ;
        tic ;
        [ t3 , COE3 ] = VoP( tspan , r0v0 , mu , tol , forces3 , A , m , epoch ) ;
        tmsise = toc ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        fprintf( 'The time to complete the exponential model of drag takes %f seconds for 120 days \n' , texp )
        fprintf( 'The time to complete the NRLMSISE-00 model of drag takes %f seconds for 120 days \n' , tmsise )

   
        figure
        hold on 
        plot( t2*s2d , COE2( :,2 )*r2d )
        plot( t3*s2d , COE3( :,2 )*r2d )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        legend( 'Exponential' , 'NRLMSISE' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,3 ) )
        plot( t3*s2d , COE3( :,3 ) )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        legend( 'Exponential' , 'NRLMSISE' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,4 )*r2d )
        plot( t3*s2d , COE3( :,4 )*r2d )
        hold off
        title( 'RAAN' )
        legend( 'Exponential' , 'NRLMSISE' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,5 )*r2d )
        plot( t3*s2d , COE3( :,5 )*r2d )
        hold off
        title( 'Argument of Perigee [degrees]' )
        legend( 'Exponential' , 'NRLMSISE' )
        
        figure
        hold on 
        plot( t2*s2d , COE2( :,8 )-re , 'b' )
        plot( t3*s2d , COE3( :,8 )-re ,'r' )
        plot( t2*s2d , COE2( :,9 )-re , 'b' )
        plot( t3*s2d , COE3( :,9 )-re , 'r')
        hold off
        title( 'Altitude [km]' )
        legend( 'Exponential' , 'NRLMSISE' )
        
        check = 'The NRLMSISE-00 model shows the satellite deorbiting faster than the \nexponential model \n' ;
        fprintf( check )
        
    end

%% Functions
function [ t , tadj , rv , COES ] = TwoBody( tspan , r0v0 , mu , tol , sf )
fprintf( "Starting Two Body \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol ) ;
    [ t , rv ] = ode45( @TwoBodyForceFun , tspan , r0v0 , opts , mu ) ;
    ii = 2 ;
    tadj(1:2) = t(1:2) ;
    COES(1,:) = state2COE( rv(1,:) , mu ) ;
    while tadj(ii-1) < t( floor( end/2 - sf(1) ) ) 
        index1 = (ii-1)*sf(1) ;
        COES(ii,:) = state2COE( rv(index1,:) , mu ) ;
        tadj(ii) = t( index1 ) ;
        ii = ii + 1 ;
    end
    jj = 0 ;
    while tadj(ii + jj-1) < t( floor( end*3/4 - sf(2) ) ) 
        index2 = index1 + (jj)*sf(2);
        COES(ii+jj,:) = state2COE( rv( index2,:) , mu ) ;
        tadj(ii+jj) = t( index2 ) ;
        jj = jj + 1 ;
    end
    jj = jj + ii ;
    kk = 0 ;
    while tadj( jj + kk-1 ) < t( floor( end - sf(3) ) ) 
        index3 = index2 + (kk)*sf(3);
        COES(jj+kk,:) = state2COE( rv( index3,:) , mu ) ;
        tadj(jj+kk) = t( index3 ) ;
        kk = kk + 1 ;
    end
    
    
    
    function drdv = TwoBodyForceFun( t , rv , mu )
        r = rv(1:3) ;
        v = rv(4:6) ;

        dr = v ;
        dv = ( -mu / norm(r)^3 )*r ;
        drdv = [ dr ; dv ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , rv , mu )
        r = norm( rv(1:3) ) ;
        value = r - ( 6378 + 100 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end

end

function [ t , tadj , rv , COES] = Cowell( tspan , r0v0 , mu , tol , forces , A , m , sf )
fprintf( "Starting Cowell \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    [ t , rv ] = ode45( @CowellForceFun , tspan , r0v0 , opts , mu , forces , A , m ) ;
    ii = 2 ;
    tadj(1:2) = t(1:2) ;
    COES(1,:) = state2COE( rv(1,:) , mu ) ;
    while tadj(ii-1) < t( floor( end/2 - sf(1) ) ) 
        index1 = (ii-1)*sf(1) ;
        COES(ii,:) = state2COE( rv(index1,:) , mu ) ;
        tadj(ii) = t( index1 ) ;
        ii = ii + 1 ;
    end
    jj = 0 ;
    while tadj(ii + jj-1) < t( floor( end*3/4 - sf(2) ) ) 
        index2 = index1 + (jj)*sf(2);
        COES(ii+jj,:) = state2COE( rv( index2,:) , mu ) ;
        tadj(ii+jj) = t( index2 ) ;
        jj = jj + 1 ;
    end
    jj = jj + ii ;
    kk = 0 ;
    while tadj( jj + kk-1 ) < t( floor( end - sf(3) ) ) 
        index3 = index2 + (kk)*sf(3);
        COES(jj+kk,:) = state2COE( rv( index3,:) , mu ) ;
        tadj(jj+kk) = t( index3 ) ;
        kk = kk + 1 ;
    end
% tadj = 0 ;
% COES = 0 ;
    
    function drdv = CowellForceFun( t , rv , mu , forces , A , m )

        r = rv(1:3) ;
        v = rv(4:6) ;
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
            if contains( forces , "nbody" )
                apn = zeros(3,1) ; 
            else
                apn = zeros(3,1) ;
            end
            if contains( forces , "srp" )
                aps = zeros(3,1) ; 
            else
                aps = zeros(3,1) ;
            end
            ap = apd + aph + apn + aps ;
        dr = v ;
        dv = ( -mu / norm(r)^3 )*r + ap ;
        drdv = [ dr ; dv ] ;
    end

    function [ value , isterminal , direction ] = Reentry( t , rv , mu , forces , A , m )
        r = norm( rv(1:3) ) ;
        value = r - ( 6378 + 100 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end
end

% function [ t , tadj , r , v , COES ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m , sf)
% fprintf( "Starting Encke \n" )
%     dr = zeros(3,1) ;
%     eps = 0 ;
%     f = 0 ;
%     rp = r0 ;
%     vp = v0 ;
%     r = r0 ;
%     v = v0 ;
%     t = tspan(1) ;
%     ii = 2 ;
%     da = zeros(3,1) ;
%     dv = zeros(3,1) ;
%     while t(ii-1) < tspan(2) && norm(r(:,ii-1)) >= (6378 + 100 ) 
% %         
%         [ r(:,ii) , v(:,ii) ] = kepler( r(:,ii-1) , v(:,ii-1) , dt ) ; %from Vallado
%         eps = dot( r(:,ii) , dr )/norm( r(:,ii) )^2 ;
%         if dr ~= zeros(3,1)
%             f = ( 1/eps )*( 1 - ( 1/ ( 1 - 2*eps )^1.5 ) ) ;
%         else
%             f = 0 ;
%         end
%             if contains( forces , "drag" )
%                 apd = DragAcceleration( r(:,ii) , v(:,ii) , A , m ) ;
%             else
%                 apd = zeros(3,1) ;
%             end
%             if contains( forces , "gravity" )
%                 if contains( forces , "J2" )
%                     aph = J2accel( r ) ; 
%                 end
%                 if contains( forces , "J3" )
%                     aph = aph + J3accel( r ) ; 
%                 end
%             else 
%                 aph = zeros(3,1) ; 
%             end
%             if contains( forces , "nbody" )
%                 apn = zeros(3,1) ; 
%             else
%                 apn = zeros(3,1) ;
%             end
%             if contains( forces , "srp" )
%                 aps = zeros(3,1) ; 
%             else
%                 aps = zeros(3,1) ;
%             end
%             ap = apd + aph + apn + aps ;
%         da = ap + ( mu/norm(r(:,ii))^3 )*( f*eps*rp - dr ) ;
%         dv = da*dt + dv ;
%         dr = .5*da*dt^2 + dv*dt + dr ;
% %         if abs(norm(dr)/norm(rp)) > 1e-2 
%             rp = r(:,ii) + dr ;
%             vp = v(:,ii) + dv ;
%             r(:,ii) = rp ;
%             v(:,ii) = vp ;
%             dr = zeros(3,1) ;
%             dv = zeros(3,1) ;
% %         else 
% %             rp = r(:,ii) + dr ;
% %             vp = v(:,ii) + dv ;
% %         end
%         t(ii) = t(ii-1) + dt ;
%         ii = ii + 1 ;
%     end
% %     tadj = 0 ;
% % COES = 0 ;
%     
%     
%     
%     ii = 2 ;
%     tadj(1:2) = t(1:2) ;
%     COES(1,:) = state2COE( [ r(:,1) ; v(:,1) ] , mu ) ;
%     while tadj(ii-1) < t( floor( end/2 - sf(1) ) ) 
%         index1 = (ii-1)*sf(1) ;
%         COES(ii,:) = state2COE( [ r(:,index1) ; v(:,index1) ] , mu ) ;
%         tadj(ii) = t( index1 ) ;
%         ii = ii + 1 ;
%     end
%     jj = 0 ;
%     while tadj(ii + jj-1) < t( floor( end*3/4 - sf(2) ) ) 
%         index2 = index1 + (jj)*sf(2);
%         COES(ii+jj,:) = state2COE( [ r(:,index2) ; v(:,index2) ] , mu ) ;
%         tadj(ii+jj) = t( index2 ) ;
%         jj = jj + 1 ;
%     end
%     jj = jj + ii ;
%     kk = 0 ;
%     while tadj( jj + kk-1 ) < t( floor( end - sf(3) ) ) 
%         index3 = index2 + (kk)*sf(3);
%         COES(jj+kk,:) = state2COE( [ r(:,index3) ; v(:,index3) ] , mu ) ;
%         tadj(jj+kk) = t( index3 ) ;
%         kk = kk + 1 ;
%     end
%     
% end

function [ t , COES ] = VoP( tspan , r0v0 , mu , tol , forces , A , m , epoch ) 
    COESo = state2COE( r0v0 , mu ) ;
    fprintf( "Starting VoP \n" )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol , 'Events' , @Reentry ) ;
    [ t , COES ] = ode45( @VoPpropagation , tspan , COESo , opts , mu , forces , A , m , epoch) ;
    COES(:,7) = ( COES(:,1).^2 ./ mu ).*( 1 ./ ( 1 - COES(:,3).^2 ) ) ;
    COES(:,8) = COES(:,7).*(1-COES(:,3)) ;
    COES(:,9) = COES(:,7).*(1+COES(:,3)) ;
    function dCOES = VoPpropagation( t , COES , mu , forces , A , m , epoch )
        
        
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
                if contains( forces , "NRLMSISE" )
                    [ epochnew , yds ] = epochUpdate( epoch , t ) ;
                    apd = DragNRLMSISE( r , v , A , m , epochnew , yds ) ;
                else
                    apd = DragAcceleration( r , v , A , m ) ;
                end
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
            if contains( forces , "nbody" )
                apn = zeros(3,1) ; 
            else
                apn = zeros(3,1) ;
            end
            if contains( forces , "srp" )
                aps = zeros(3,1) ; 
            else
                aps = zeros(3,1) ;
            end
            apeci = apd + aph + apn + aps ;
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

    function [ value , isterminal , direction ] = Reentry( t , COES , mu , forces , A , m , epoch )
        a = ( COES(1).^2 ./ mu ).*( 1 ./ ( 1 - COES(3).^2 ) ) ;
        rp = a.*(1-COES(3)) ;
        value = rp - ( 6378 + 100 ) ;
        isterminal = 1 ;
        direction = 0 ;
    end
end

function aDrag = DragNRLMSISE( r , v , A , m , epoch , yds)
    lla = eci2lla( r'*1e3 , epoch ) ;
    [ ~ , rho ] = atmosnrlmsise00( lla(3) , lla(1) , lla(2) , yds(1) , yds(2) , yds(3) ) ; 
    CD = 2.2 ;
    wearth = [ 0 ; 0 ; 7.29211e-5 ] ;
    vrel = v - cross( wearth , r ) ;
    aDrag = ( -.5*((CD*A)/m)*rho(6)*norm(vrel)*vrel*1e6 )*1e-3 ;
end

function aDrag = DragAcceleration( r , v , A , m )
    CD = 2.2 ;
    wearth = [ 0 ; 0 ; 7.29211e-5 ] ;
        h0 = [ 0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000 ]' ;
        rho0 = [ 1.225 , 3.899e-2 ,  1.774e-2 , 3.972e-3 , 1.057e-3 , 3.206e-4 , 8.770e-5 , 1.905e-5 , 3.396e-6 , 5.297e-7 , 9.661e-8 , ...
        2.438e-8 , 8.484e-9 , 3.845e-9 , 2.070e-9 , 5.464e-10 , 2.789e-10 , 7.248e-11 , 2.418e-11 , 9.518e-12 , 3.725e-12 , ...
        1.585e-12 , 6.967e-13 , 1.454e-13 , 3.614e-14 , 1.170e-14 , 5.245e-15 , 3.019e-15 ]' ;
        H = [ 7.249 , 6.349 , 6.682 , 7.554 , 8.382 , 7.714 , 6.549 , 5.799 , 5.382 , 5.877 , 7.263 , 9.473 , 12.636 , 16.149 , ... 
        22.523 , 29.740 , 37.105 , 45.546 , 53.628 , 53.298 , 58.515 , 60.828 , 63.822 , 71.835 , 88.667 , 124.64 , 181.05 , 286.00 ]' ;
    hellp = norm(r) - 6378 ;
    index = max( find( h0 < hellp ) ) ;
    rho = rho0(index)*exp( ( - ( hellp - h0(index)) )/ H(index) ) ;
    vrel = v - cross( wearth , r ) ;
    aDrag = ( -.5*((CD*A)/m)*rho*norm(vrel)*vrel*1e6 )*1e-3 ;
end

function COES = state2COE( state, mu )
% all angles output in degrees
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
    Kh = [ 0 0 1 ] ; % K hat
    R = state(1:3) ;
    V = state(4:6) ;
    r = norm( R ) ;
    v = norm( V ) ;
    vr = dot( R , V )/r ; % radial velocity
    H = cross( R , V ) ; % specific angular momentum
    h = norm( H ) ; % specific angular momentum
    inc = acos( H(3)/h ) ; %inclination
    ECC = (1/mu)*( ( v^2 - mu/r )*R - r*vr*V ) ; %eccentricity vector
    ecc = norm( ECC ) ; % eccentricity
    N = cross( Kh , H ) ; % Node line
    n = norm( N ) ;
if n ~= 0
    RAAN = acos(N(1)/n) ; %Right ascension of ascending node
    if N(2) < 0 
        RAAN = 2*pi - RAAN ; %Right ascension of ascending node
    end
else
    RAAN = 0 ;
end 
    
if n ~= 0 
    if ecc >= 0 
        omega = acos(dot(N,ECC)/(n*ecc)) ; % Argument of perigee
        if ECC(3) < 0 
            omega = 2*pi - omega ; % Argument of perigee
        end
    else
        omega = 0 ;
    end
else
    omega = 0 ;
end
    
if ecc > 0
    theta = acos( dot( ECC , R )/( ecc*r ) ) ;         
    if vr < 0 
        theta = 2*pi - theta ; 
    end
else
    cp = cross( N , R ) ;
    if cp(3) >= 0
        theta = acos( dot( N , R )/( n*r ) ) ;
    else
        theta = 2*pi - acos( dot( N , R )/( n*r ) ) ;
    end
end

    a = (h^2)/( mu*( 1 - ecc^2 ) ) ; % semi-major axis
    rp = a*( 1 - ecc ) ;
    ra = a*( 1 + ecc ) ;

    
    COES = [ h , inc , ecc , RAAN , omega , theta , a , rp , ra ] ;
end

function [ r , v ] = coes2state( COES , mu )
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
h = COES(1) ;
inc = COES(2) ;
ecc = COES(3) ;
RAAN = COES(4) ;
omega = COES(5) ;
theta = COES(6) ;

    r_peri = (h^2/mu) * ( 1/( 1 + ecc*cos(theta) ) ) * [ cos( theta ) ; sin( theta ) ; 0 ] ;
    v_peri = (mu/h) * [ -sin( theta ) ; ecc+cos(theta) ; 0 ] ;

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

function aP = J2accel( r )
mu = 398600 ;
J2 = -1.08262617385222e-3 ;
Re = 6378 ;
aP = zeros(3,1) ;
rmag = norm(r) ;
    aP(1) = -3*J2*mu*Re^2*r(1)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
    aP(2) = -3*J2*mu*Re^2*r(2)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
    aP(3) = -3*J2*mu*Re^2*r(3)/(2*rmag^5)*(3-5*r(3)^2/rmag^2) ;
end

function aP = J3accel( r )
mu = 398600 ;
J3 = 2.53241051856772e-6 ;
Re = 6378 ;
aP = zeros(3,1) ;
rmag = norm(r) ;
    aP(1) = -5*J3*mu*Re^3*r(1)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
    aP(2) = -5*J3*mu*Re^3*r(2)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
    aP(3) = -5*J3*mu*Re^3/(2*rmag^7)*(6*r(3)^2-7*r(3)^4/rmag^2 - (3/5)*rmag^2) ;
end
% 
%     function [r, v] =  kepler  ( ro, vo, dtseco )
%           % ------------------------------------------------------------------------------
%     %
%     %                           function kepler
%     %
%     %  this function solves keplers problem for orbit determination and returns a
%     %    future geocentric equatorial (ijk) position and velocity vector.  the
%     %    solution uses universal variables.
%     %
%     %  author        : david vallado                  719-573-2600   22 jun 2002
%     %
%     %  revisions
%     %    vallado     - fix some mistakes                             13 apr 2004
%     %
%     %  inputs          description                    range / units
%     %    ro          - ijk position vector - initial  km
%     %    vo          - ijk velocity vector - initial  km / s
%     %    dtsec       - length of time to propagate    s
%     %
%     %  outputs       :
%     %    r           - ijk position vector            km
%     %    v           - ijk velocity vector            km / s
%     %    error       - error flag                     'ok', ...
%     %
%     %  locals        :
%     %    f           - f expression
%     %    g           - g expression
%     %    fdot        - f dot expression
%     %    gdot        - g dot expression
%     %    xold        - old universal variable x
%     %    xoldsqrd    - xold squared
%     %    xnew        - new universal variable x
%     %    xnewsqrd    - xnew squared
%     %    znew        - new value of z
%     %    c2new       - c2(psi) function
%     %    c3new       - c3(psi) function
%     %    dtsec       - change in time                 s
%     %    timenew     - new time                       s
%     %    rdotv       - result of ro dot vo
%     %    a           - semi or axis                   km
%     %    alpha       - reciprocol  1/a
%     %    sme         - specific mech energy           km2 / s2
%     %    period      - time period for satellite      s
%     %    s           - variable for parabolic case
%     %    w           - variable for parabolic case
%     %    h           - angular momentum vector
%     %    temp        - temporary real*8 value
%     %    i           - index
%     %
%     %  coupling      :
%     %    mag         - magnitude of a vector
%     %    findc2c3    - find c2 and c3 functions
%     %
%     %  references    :
%     %    vallado       2004, 95-103, alg 8, ex 2-4
%     %
%     % [r, v] =  kepler  ( ro, vo, dtsec );
%     % ------------------------------------------------------------------------------
%         
%     %function [r,v,errork] =  kepler  ( ro,vo, dtseco, fid );
% 
%     % -------------------------  implementation   -----------------
%     % set constants and intermediate printouts
% mu = 398600.4418 ;
%     numiter    =    50;
% small = 1e-10 ;
% twopi = 2*pi ;
%     
% 
%     % --------------------  initialize values   -------------------
%     znew  = 0.0;
%     dtsec = dtseco;
% 
%     if ( abs( dtseco ) > small )
%         magro = mag( ro );
%         magvo = mag( vo );
%         rdotv= dot( ro,vo );
% 
%         % -------------  find sme, alpha, and a  ------------------
%         sme= ( (magvo^2)*0.5  ) - ( mu /magro );
%         alpha= -sme*2.0/mu;
% 
%         if ( abs( sme ) > small )
%             a= -mu / ( 2.0 *sme );
%         else
%             a= infinite;
%         end
%         if ( abs( alpha ) < small )   % parabola
%             alpha= 0.0;
%         end
% 
% 
%         % ------------   setup initial guess for x  ---------------
%         % -----------------  circle and ellipse -------------------
%         if ( alpha >= small )
%             period= twopi * sqrt( abs(a)^3.0/mu  );
%             % ------- next if needed for 2body multi-rev ----------
%             if ( abs( dtseco ) > abs( period ) )
%                 % including the truncation will produce vertical lines that are parallel
%                 % (plotting chi vs time)
%                 dtsec = rem( dtseco, period );
%             end
%             xold = sqrt(mu)*dtsec * alpha;
%         else
%             % --------------------  parabola  ---------------------
%             if ( abs( alpha ) < small )
%                 h = cross( ro,vo );
%                 magh = mag(h);
%                 p= magh*magh/mu;
%                 s= 0.5  * (halfpi - atan( 3.0 *sqrt( mu / (p*p*p) )* dtsec ) );
%                 w= atan( tan( s )^(1.0 /3.0 ) );
%                 xold = sqrt(p) * ( 2.0 *cot(2.0 *w) );
%                 alpha= 0.0;
%             else
%                 % ------------------  hyperbola  ------------------
%                 temp= -2.0 * mu * dtsec / ...
%                     ( a*( rdotv + sign(dtsec)*sqrt(-mu*a)* ...
%                     (1.0 -magro*alpha) ) );
%                 xold= sign(dtsec) * sqrt(-a) *log(temp);
%             end
%         end
% 
%         ktr= 1;
%         dtnew = -10.0;
%         % conv for dtsec to x units
%         tmp = 1.0 / sqrt(mu);
% 
%         while ((abs(dtnew*tmp - dtsec) >= small) && (ktr < numiter))
%             xoldsqrd = xold*xold;
%             znew     = xoldsqrd * alpha;
% 
%             % ------------- find c2 and c3 functions --------------
%             [c2new, c3new] = findc2c3( znew );
% 
%             % ------- use a newton iteration for new values -------
%             rval = xoldsqrd*c2new + rdotv*tmp *xold*(1.0 -znew*c3new) + ...
%                 magro*( 1.0  - znew*c2new );
%             dtnew= xoldsqrd*xold*c3new + rdotv*tmp*xoldsqrd*c2new + ...
%                 magro*xold*( 1.0  - znew*c3new );
% 
%             % ------------- calculate new value for x -------------
%             temp1 = ( dtsec*sqrt(mu) - dtnew ) / rval;
%             xnew = xold + temp1;
% 
%             % ----- check if the univ param goes negative. if so, use bissection
%             if xnew < 0.0
%                 xnew = xold*0.5;
%             end
%   
%             ktr = ktr + 1;
%             xold = xnew;
%         end
% 
%         if ( ktr >= numiter )
%             errork= 'knotconv';
%             fprintf(1,'kep not conv in %2i iter %11.3f \n',numiter, dtseco );
%             for i= 1 : 3
%                 v(i)= 0.0;
%                 r(i)= v(i);
%             end
%         else
%             % --- find position and velocity vectors at new time --
%             xnewsqrd = xnew*xnew;
%             f = 1.0  - ( xnewsqrd*c2new / magro );
%             g = dtsec - xnewsqrd*xnew*c3new/sqrt(mu);
% 
%             for i= 1 : 3
%                 r(i)= f*ro(i) + g*vo(i);
%             end
%             magr = mag( r );
%             gdot = 1.0  - ( xnewsqrd*c2new / magr );
%             fdot = ( sqrt(mu)*xnew / ( magro*magr ) ) * ( znew*c3new-1.0  );
%             for i= 1 : 3
%                 v(i)= fdot*ro(i) + gdot*vo(i);
%             end
% 
%             temp= f*gdot - fdot*g;
%             if ( abs(temp-1.0 ) > 0.00001  )
%                 errork= 'fandg';
%             end
% 
%         end  % if
%     else
%         % ----------- set vectors to incoming since 0 time --------
%         for i=1:3
%             r(i)= ro(i);
%             v(i)= vo(i);
%         end
%     end
%     
%     end
% 

function [ t , tadj , r , v , COES ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m , sf)
    dr = zeros(3,1) ;
    eps = 0 ;
    f = 0 ;
    rp = r0 ;
    vp = v0 ;
    r = r0 ;
    v = v0 ;
    t = tspan(1) ;
    ii = 2 ;
    da = zeros(3,1) ;
    dv = zeros(3,1) ;
    while t(ii-1) < tspan(2) && norm(r(:,ii-1)) >= (6378 + 100 ) 
%         
        [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt ) ; %from Vallado
        eps = dot( r(:,ii) , dr )/norm( r(:,ii) )^2 ;
        if dr ~= zeros(3,1)
            f = ( 1/eps )*( 1 - ( 1/ ( 1 - 2*eps )^1.5 ) ) ;
        else
            f = 0 ;
        end
            if contains( forces , "drag" )
                apd = DragAcceleration( r(:,ii) , v(:,ii) , A , m ) ;
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
            if contains( forces , "nbody" )
                apn = zeros(3,1) ; 
            else
                apn = zeros(3,1) ;
            end
            if contains( forces , "srp" )
                aps = zeros(3,1) ; 
            else
                aps = zeros(3,1) ;
            end
            ap = apd + aph + apn + aps ;
        da = ap + ( mu/norm(r(:,ii))^3 )*( f*eps*rp - dr ) ;
        dv = da*dt + dv ;
        dr = .5*da*dt^2 + dv*dt + dr ;
%         if abs(norm(dr)/norm(rp)) > 1e-2 
            rp = r(:,ii) + dr ;
            vp = v(:,ii) + dv ;
            r(:,ii) = rp ;
            v(:,ii) = vp ;
            dr = zeros(3,1) ;
            dv = zeros(3,1) ;
%         else 
%             rp = r(:,ii) + dr ;
%             vp = v(:,ii) + dv ;
%         end
        t(ii) = t(ii-1) + dt ;
        ii = ii + 1 ;
    end

    tadj = [ t(1) , t(2) ] ;
    COES(1,:) = state2COE( [ r(:,1) , v(:,1) ] , mu ) ;
    sf(1) = floor(length(t)/1e4) ;
    sf(2) = floor(length(t)/1e4) ;
    sf(3) = floor(length(t)/1e4) ;
    while tadj(ii-1) < t( floor( end/2 - sf(1) ) ) 
        index1 = (ii-1)*sf(1) ;
        COES(ii,:) = state2COE( rv(index1,:) , mu ) ;
        tadj(ii) = t( index1 ) ;
        ii = ii + 1 ;
    end
    jj = 0 ;
    while tadj(ii + jj-1) < t( floor( end*3/4 - sf(2) ) ) 
        index2 = index1 + (jj)*sf(2);
        COES(ii+jj,:) = state2COE( rv( index2,:) , mu ) ;
        tadj(ii+jj) = t( index2 ) ;
        jj = jj + 1 ;
    end
    jj = jj + ii ;
    kk = 0 ;
    while tadj( jj + kk-1 ) < t( floor( end - sf(3) ) ) 
        index3 = index2 + (kk)*sf(3);
        COES(jj+kk,:) = state2COE( rv( index3,:) , mu ) ;
        tadj(jj+kk) = t( index3 ) ;
        kk = kk + 1 ;
    end
    
    
end






end