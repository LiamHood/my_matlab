function [ r , v , RMS , P , AtWA , AtWb ] = DCODrazel2rv( obs , obserr , lat , long , alt , UTC )
% Performs differential correction orbit determination using Weighted Least
% Squares technique. The observations for this function should be razel
% columns. 
    mu = 398600 ;
    d2s = 86400 ;
    rho0 = obs{1} ;
    az0 = obs{2} ;
    el0 = obs{3} ;
    rhoerr = obserr(1) ;
    azerr = obserr(2) ;
    elerr = obserr(3) ;
    n = length(el0) ; 
    
    
    W = zeros( 3 ) ;
    W(1,1) = 1/( rhoerr^2 ) ;
    W(2,2) = 1/( azerr^2 ) ;
    W(3,3) = 1/( elerr^2 ) ;
%     W = W/W(1,1) ;
%     W = W/norm(W) ;
    
    RMS0 = 1 ;
        for ii = 1:n %find position vector from observation
            [ r(:,ii) ] = razel2r( rho0(ii) , az0(ii) , el0(ii) , lat , long , alt , UTC(ii,:) ) ;
            JD(ii) = juliandate( UTC(ii,:) ) ;
            obso(:,ii) = [ rho0(ii) ; az0(ii) ; el0(ii) ] ;
        end
        for ii = 2:(n-1) % find velocities for all but first and last observation
            v(:,ii) = Gibbs( r(:,ii-1) , r(:,ii) , r(:,ii+1) , JD(ii-1) , JD(ii) , JD(ii+1) , mu ) ;
        end
        for ii = 2:n-1 % put all states back to initial epoch
            tspan = [ 0 , JD(1) - JD(ii)]*d2s ;
            [ ~ , rback , vback ] = TwoBody( tspan , r(:,ii) , v(:,ii) , mu , 1e-8 ) ;
            r0(:,ii-1) = rback(:,end) ;
            v0(:,ii-1) = vback(:,end) ;
        end
        r0avg = zeros( 3 , 1 ) ;
        v0avg = zeros( 3 , 1 ) ;
        for ii = 1:3 % find average state at initial epoch
            r0avg(ii) = mean( r0(ii,:) ) ;
            v0avg(ii) = mean( v0(ii,:) ) ;
        end
        rnom = r0avg ;
        vnom = v0avg ;
        xnom0 = [ rnom ; vnom ] ;
        xnom0 = [5975.29040000000;2568.64000000000;3120.58450000000;-3.98384600000000;-2.07115900000000;-5.91709500000000] ;
        xnom = zeros( 6 , n ) ;
        
    tol = 1e-3 ;
    err = 1 ;
    while err >= tol % run until RMS changes by less that %0.1
%         A = zeros( 3 , 6 ) ;
%         btil = zeros( n , 1 ) ;
        AtWA = zeros( 6 , 6 ) ;
        AtWb = zeros( 6 , 1 ) ;
        for ii = 1:n %loop for all observations
            tspan = [ 0 , JD(ii) - JD(1) ]*d2s ;
            if tspan(1) ~= tspan(2)
                [ ~ , rnomi , vnomi ] = TwoBody( tspan , xnom0(1:3) , xnom0(4:6) , mu , 1e-8 ) ;
                xnom(:,ii) = [ rnomi(:,end) ; vnomi(:,end) ] ;
            else
                xnom(:,1) = xnom0 ;
                rnomi = xnom0(1:3) ;
                vnomi = xnom0(4:6) ;
            end
            [ rhoni , azni , elni , ~ , ~ , ~ , ~ ] = RAZEL( rnomi(:,end) , vnomi(:,end) , UTC(ii,:) , lat , long , alt ) ;
            btil = [ rho0(ii) - rhoni ; az0(ii) - azni ; el0(ii) - elni ]  ;
            for jj = 1:6 % finite differencing for every element
                xmod = xnom0(:) ;
                delement = xnom0(jj)*.0001 ;
                xmod(jj) = xnom0(jj) + delement ;
                if tspan(2) ~= 0 
                    [ ~ , rmodi , vmodi ] = TwoBody( tspan , xmod(1:3) , xmod(4:6) , mu , 1e-8 ) ;
                else
                    rmodi = xmod(1:3) ;
                    vmodi = xmod(4:6) ;
                end
                [ rhomi , azmi , elmi , ~ , ~ , ~ , ~ ] = RAZEL( rmodi(:,end) , vmodi(:,end) , UTC(ii,:) , lat , long , alt ) ;
                A(:,jj) = [ ( rhomi - rhoni )/delement ; ( azmi - azni )/delement ; ( elmi - elni )/delement ] ;
            end 
            AtWAi = A'*W*A ;
            AtWA = AtWA + AtWAi ;
            AtWbi = A'*W*btil ;
            AtWb = AtWb + AtWbi ;
        end
        P = pinv(AtWA) ;
        delx = P*AtWb ;
        RMS = sqrt( ( btil'*W*btil )/( 3*n ) ) ;
        err = abs( RMS - RMS0 )/RMS ;
        if err >= tol 
            xnom0 = xnom0 + delx ;
        end
        RMS0 = RMS ;
    end
    r = xnom0(1:3) ;
    v = xnom0(4:6) ;
end