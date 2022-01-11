function [ r , v , P ] = SBF( newobs , obserr , lat , long , alt , UTC , UTC0 , rold , vold , AtWAold , AtWbold )
     mu = 398600 ;
    d2s = 86400 ;
    rho0 = newobs{1} ;
    az0 = newobs{2} ;
    el0 = newobs{3} ;
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
    
    xnom0 = [ rold ; vold ] ;
    JD0 = juliandate( UTC0 ) ;
    
    
    AtWA = zeros( 6 , 6 ) ;
        AtWb = zeros( 6 , 1 ) ;
        for ii = 1:n %loop for all observations
            JD(ii) = juliandate( UTC(ii,:) ) ;
            tspan = [ 0 , JD(ii) - JD0 ]*d2s ;
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
        
        delx = pinv( AtWA + AtWAold )*( AtWb + AtWbold ) ;
        P = pinv( AtWA + AtWAold ) ;
        r = rold + delx(1:3) ;
        v = vold + delx(4:6) ;
end