function [t, state, COES] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m, tol, Re)
    
    dr = zeros(3,1) ;
    eps = 0;
    f = 0;
    rp = r0 ;
    vp = v0;
    r = r0 ;
    v = v0 ;
    COES = state2coes([r; v], mu );
    t = tspan(1) ;
    ii = 2 ;
    dv = zeros(3,1) ;
    while t(ii-1) < tspan(2) && COES(ii-1, 8) >= (Re + 120 ) 
        t(end)/(3600*24*365)
        % propagate unperturbed orbit
        [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt , mu  ) ; 
        
        % find eps and f for encke
        eps = dot( r(:,ii) , dr )/norm( r(:,ii) )^2 ;
        if eps ~= 0
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
    
        % Find COES so that the function developed for VOP work
        COES(ii, :) = state2coes([r(:,ii-1); v(:,ii-1)], mu ); % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        
        ap = apd + aph;
        da = ap + ( mu/norm(r(:,ii))^3 )*( f*eps*rp - dr ) ; % difference in acceleration
        dv = da*dt + dv ; % diff velocity
        dr = .5*da*dt^2 + dv*dt + dr ; % diff position
        rp = r(:,ii) + dr ; % position perturbed
        vp = v(:,ii) + dv ; % velocity perturbed
%         track_dr(ii) = norm(dr);        
        if norm(dr)/norm(rp) > tol
            r(:,ii) = rp ;
            v(:,ii) = vp ;
            dr = zeros(3,1);
            dv = zeros(3,1);
        end  
        
        t(ii) = t(ii-1) + dt ;
        ii = ii + 1 ;
    end
    state = [r; v];
end