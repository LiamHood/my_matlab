function [ r , v , RMS , P ] = WLLSrazel2rv( razel , razelerr , lla , UTC )
    rho0 = razel{1} ;
    az0 = razel{2} ;
    el0 = razel{3} ;
    rhoerr = razelerr{1} ;
    azerr = razelerr{2} ;
    elerr = razelerr{3} ;
    lat = lla(1) ;
    long = lla(2) ;
    alt = lla(3) ;
    
    mu = 398600 ;
     W = zeros( n ) ;
     W(1,1) = 1/( rhoerr^2 ) ;
     W(2,2) = 1/( azerr^2 ) ;
     W(3,3) = 1/( elerr^2 ) ;
     W = W/norm(W) ;
            
        n = length(el0) ;
        for ii = 1:n
            [ r(:,ii) ] = razel2r( rho0(ii) , az0(ii) , el0(ii) , lat , long , alt , UTC(ii,:) ) ;
            JD(ii) = juliandate( UTC(ii,:) ) ;
            obso(:,ii) = [ rho0(ii) ; az0(ii) ; el0(ii) ] ;
        end
        for ii = 2:(n-1)
            v(:,ii) = Gibbs( r(:,ii-1) , r(:,ii) , r(:,ii+1) , JD(ii-1) , JD(ii) , JD(ii+1) , mu ) ;
        end
        dt = ( JD(n) - JD(1) )*86400 ;
        [ v(:,1) , v(:,n) ] = Lambert_Battin( r(:,1) , r(:,n) , dt , 1 , mu ) ;
        r0(:,1) = r(:,1) ;
        v0(:,1) = r(:,1) ;
        for ii = 2:n
            tspan = [ 0 , JD(1) - JD(ii)]*86400 ;
            [ ~ , rback , vback ] = TwoBody( tspan , r(:,ii) , v(:,ii) , mu , 1e-8 ) ;
            r0(:,ii) = rback(:,end) ;
            v0(:,ii) = vback(:,end) ;
        end
        r0avg = zeros( 3 , 1 ) ;
        v0avg = zeros( 3 , 1 ) ;
        for ii = 1:3
            r0avg(ii) = mean( r0(ii,:) ) ;
            v0avg(ii) = mean( v0(ii,:) ) ;
        end
        rnom = r0avg ;
        vnom = v0avg ;
        while err >= tol
            [ rhoc , azc , elc , ~ , ~ , ~ , ~ ] = RAZEL( rnom , vnom , UTC(1,:) , lat , long , alt ) ;
            obsci = [ rhoc ; azc ; elc ] ;
            obsnom(:,1) = obsci ;
            xnom(:,1) = [ rnom ; vnom ] ;
            for ii= 2:n
                tspan = [ 0 , JD(ii) - JD(1) ]*86400 ;
                [ ~ , rnom , vnom ] = TwoBody( tspan , xnom(1:3,1) , xnom(4:6,1) , mu , 1e-8 ) ;
                xnom(:,ii) = [ rnom ; vnom ] ;
                [ rhoc , azc , elc , ~ , ~ , ~ , ~ ] = RAZEL( rnom , vnom , UTC(ii,:) , lat , long , alt ) ;
                obsci = [ rhoc ; azc ; elc ] ;
                obsnom(:,ii) = obsci ;
            end
            
            for jj = 1:length( xnom(1,:) )
                xmod(:,1) = xnom(:,1) ;
                del = xnom(jj,1)*(.001) ;
                xmod(jj,1) = xnom(jj,1) + del ;
                rmod = xmod(1:3,1) ;
                vmod = xmod(4:6,1) ;
                [ rhom , azm , elm , ~ , ~ , ~ , ~ ] = RAZEL( rmod , vmod , UTC(1,:) , lat , long , alt ) ;
                obsmi = [ rhom ; azm ; elm ] ;
                obsm(:,1) = obsmi ;
                for ii= 2:n
                    tspan = [ 0 , JD(ii) - JD(1) ]*86400 ;
                    [ ~ , rmod , vmod ] = TwoBody( tspan , xmod(1:3,1) , xmod(4:6,1) , mu , 1e-8 ) ;
                    xmod(:,ii) = [ rmod ; vmod ] ;
                    [ rhom , azm , elm , ~ , ~ , ~ , ~ ] = RAZEL( rmod , vmod , UTC(ii,:) , lat , long , alt ) ;
                    obsmi = [ rhom ; azm ; elm ] ;
                    obsm(:,ii) = obsmi ;
                end
                res(:,jj) = obsm
                for ii = 1:n
                    A(ii,1) = 1 ;
                    A(ii,2) = 1 ;
                end
            end

           
            
        end
end