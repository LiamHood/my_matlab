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