clear ; close all ; clc ;
tol = 1e-8 ;
opts = optimoptions( 'fsolve' , 'Display' , 'iter-detailed'  , 'FunctionTolerance' , tol ,  ...
    'OptimalityTolerance' , tol , 'MaxIterations' , 1e6 , 'MaxFunctionEvaluations' , 1e6 , 'Algorithm' , 'levenberg-marquardt' ) ; % 
[ lam , F ] = fsolve( @AlfanoTransfer , 1 , opts ) ;
    [ a , R , inc , dinc , t , r , v , km2DU , s2TU , thetac , athrust , mi , dvacc , tf ] = AlfanoTransferProp(lam) ;
    rnc = r.*km2DU ;
    tdays = t.*s2TU./86400 ;
    athrustnc = athrust./(s2TU^2/km2DU) ;
    thrust = athrustnc.*1e3.*mi ;
    figure
    axis equal
    plot3( rnc(1,:) , rnc(2,:) , rnc(3,:) )
    figure
    plot( tdays , thrust )
    figure
    plot( tdays , thetac )
    
function F = AlfanoTransfer( lam )
    [ a , R , inc , dinc , ~ , ~ , ~ , ~ , ~ , ~ , ~ , ~ , ~ , ~] = AlfanoTransferProp(lam) ;
    F(1,1) = a - R ;
    F(2,1) = inc - dinc ;
end

function [ a , R , inc , dinc , t , r , v , km2DU , s2TU , thetac , athrust , mi , dvacc , tf ] = AlfanoTransferProp(lam)
    d2r = pi/180 ;
    mu = 4905 ;
    ainit = 7500 ;
    afinal = 7500 ;
    iinit = 35*d2r ;
    ifinal = 55*d2r ;
    F = 0.014 ;
    Isp = 1100 ;
    mi = 150 ;
    t = 0 ;

    R = afinal / ainit ;
    mdot = -F/(Isp*9.81) ; % F = mdot*Isp*9.81
    mds = mdot/mi ;
    dinc = ifinal - iinit ;
    athrustiNCU = F/mi*1e-3 ;
    km2DU = ainit ;
    s2TU = sqrt( ainit^3/mu ) ;
    athrust(1) = athrustiNCU*s2TU^2/km2DU ;
    a = 1 ;
    inc = 0 ;

    
    costp = cos( 0 ) ;
    steps = 1000 ;
    dt = 2*pi*sqrt(a^3)/steps ;
    ii = 2 ;
    [ r(:,1) , v(:,1) ] = coes2state( [ sqrt(a) , 0 , 1e-12 , 0 , 0 , 0 ] , 1 ) ;
    [ rf , vf ] = coes2state( [ sqrt(R) , dinc , 1e-12 , 0 , 0 , 0 ] , 1 ) ;
    hf = cross( rf , vf ) ;
    beta0 = atan3(sin(0.5 * pi * dinc), (norm(v(:,1))/norm(vf)) - cos(0.5 * pi * dinc));
        while inc < dinc
            [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt , 1  ) ; 
            h = cross( r(:,ii) , v(:,ii) ) ;
            ncurrent = cross( [ 0 ; 0 ; 1 ] , h ) ;
            if r(3,ii) > 0
                arglat = acos( dot( ncurrent , r(:,ii) )/( norm(ncurrent)*norm(r(:,ii)) ) ) ;
            else
                arglat = acos( dot( ncurrent , r(:,ii) )/( norm(ncurrent)*norm(r(:,ii)) ) ) + pi ;
            end
            athrust(ii) = athrust(ii-1)/( 1 + mds*dt ) ;
%         alpha = [ 0 , 2.467410607 , -1.907470562 , 35.892442177 , -214.672979624 , ...
%             947.773273608 , -2114.861134906 , 2271.240058672 , -1127.457440108 , ...
%             192.953875268 , 8.577733773 ] ;
%         beta = [ 1 , 0.4609698838 , 13.7756315324 , -69.1245316678 , 279.0671832500 ,...
%             -397.6628952136 , -70.0139935047 , 528.0334266841 , -324.9303836520 ...
%             20.5838245170 , 18.8165370778 ] ;
%         alphasum = 0 ;
%         betasum = 0 ;
%         x = pi/(2*lam)*sqrt( 1/a ) ;
%         z = 1/(x^2) ;
%         for jj = 1:11
%            alphasum = alphasum + alpha(jj)*z^jj ;
%            betasum = betasum + beta(jj)*z^jj ;
%         end
%         cv = alphasum*betasum ;
            cv = 1/(4*lam^2*a^2+1) ;
    %         thetap = acos( dot( r , ncommon )/norm( r ) ) ;
            sinci = costp ;
            cosci = sqrt( (1/cv) - 1 ) ;
            sinc = sinci/sqrt( sinci^2 + cosci^2 ) ;
            cosc = cosci/sqrt( sinci^2 + cosci^2 ) ;
            if isreal( sinc ) && isreal( cosc )
                thetac(ii) = atan3( sinc , cosc ) ;
%                 thetac(ii) = atan3( norm(v(:,ii)) * sin(beta0) , (norm(v(:,ii)) * cos(beta0) - norm(athrust(ii)) * dt ) ) ;
%                 beta_tmp = atan3(v1 * sin(beta0), (v1 * cos(beta0) - thracc * t));
            else
                thetac = 1 ;
            end
            if (arglat >= 0.0 && arglat < 0.5 * pi)  
                thetac(ii) = thetac(ii) ;  
            end
            if (arglat > 0.5 * pi && arglat < pi)
                thetac(ii) = thetac(ii); 
            end
            if (arglat > pi && arglat < 1.5 * pi)
               thetac(ii) = thetac(ii);
            end
            if (arglat > 1.5 * pi && arglat < 2.0 * pi) 
               thetac(ii) = thetac(ii);
            end
            antw = athrust(ii).* [ 0 ; cos( thetac(ii) ) ; sin( thetac(ii) ) ] ;
            tv = v(:,ii)/norm(v(:,ii)) ;
            wv = h/norm( h ) ;
            nv = cross( tv , wv )/( norm( tv )*norm( wv ) ) ;
            ntw2eci = [ nv , tv , wv ] ;
            aeci = ntw2eci*antw ;

                da = aeci ;
                dv = da*dt ;
                dr = .5*da*dt^2 + dv*dt ;
                rp = r(:,ii) + dr ;
                vp = v(:,ii) + dv ;
                r(:,ii) = rp ;
                v(:,ii) = vp ;
                t(ii) = t(ii-1) + dt ;
                h = cross( r(:,ii) , v(:,ii) ) ;
                ncommon = cross( hf , h ) ;
                ncommon = ncommon/norm( ncommon ) ;
                costp = dot( r(:,ii) , ncommon )/norm( r(:,ii) ) ;
                sme = norm( v(:,ii) )^2/2 - 1/norm( r(:,ii) ) ;
                a = -2*sme ;
                inc = acos( h(3) / norm( h ) ) ;
                dt = 2*pi*sqrt(a^3)/steps ;
                ii = ii + 1 ;
        end
        dvacc = 1 - sqrt( 1/R ) ;
        tf = 1/(-mds)*(1-exp(mds*dvacc/athrust(ii-1)))*s2TU ;
end
    
function ds = LongTimeScale( t , s , lam )
    a = s(1) ;
    inc = s(2) ;
    alpha = [ 0 , 2.467410607 , -1.907470562 , 35.892442177 , -214.672979624 , ...
        947.773273608 , -2114.861134906 , 2271.240058672 , -1127.457440108 , ...
        192.953875268 , 8.577733773 ] ;
    beta = [ 1 , 0.4609698838 , 13.7756315324 , -69.1245316678 , 279.0671832500 ,...
        -397.6628952136 , -70.0139935047 , 528.0334266841 , -324.9303836520 ...
        20.5838245170 , 18.8165370778 ] ;
    alphasum = 0 ;
    betasum = 0 ;
    x = pi/(2*lam)*sqrt( 1/a ) ;
    z = 1/(x^2) ;
    for ii = 1:11
       alphasum = alphasum + alpha(ii)*z^ii ;
       betasum = betasum + beta(ii)*z^ii ;
    end
    cv = alphasum*betasum ;
    [ K , E ] = ellipke( cv , 1e-12 ) ;
    ds(1,1) = (4/pi)*sqrt( a^3 )*sqrt( 1-cv )*K ;
    ds(2,1) = (2/pi)*sqrt(a)*( (1/sqrt(cv))*E + ( sqrt(cv) - 1/sqrt(cv) )*K ) ;
    
end