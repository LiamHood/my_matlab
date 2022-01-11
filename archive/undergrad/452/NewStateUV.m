% function [ r , v ] = NewStateUV( r0 , v0 , dt , mu , denomS , denomC )
% % Find new position and velocity at some time in orbit by lagrange
% % variables
% 
%     % Initial values
%     r0mag = norm( r0 ) ; % radius in km
%     v0mag = norm( v0 ) ; % speed in km/s
%     vr0 = dot( r0 , v0 ) / r0mag ; % radial speed (km/s)   
%     alpha = (2/r0mag) - ( v0mag^2 / mu ) ; % reciporical of semimajor axis
%     chi = sqrt( mu ) * dt * abs( alpha ) ; % Initial universal anomaly guess
%     z = alpha*chi(1)^2 ;    
%     [ S , C ] = StumpfCalc( z , denomS , denomC ) ;
%     
%     % Universal anomaly equations
%     
%     % Curtis
%     fun = @(chi,C,S) ( ( r0mag*vr0 ) / sqrt( mu ) )*chi^2*C + ( 1 - alpha*r0mag )*chi^3*S + r0mag*chi - sqrt( mu )*dt ;
%     fprime = @(chi,C,S) ( ( r0mag*vr0 ) / sqrt( mu ) )*chi*( 1 - alpha*chi^2*S ) + ( 1 - alpha*r0mag )*chi^2*C + r0mag ; 
%     
%         % Newtons for UV
%         ii = 1 ;
%         ratio = 1 ;
%         tol = 1e-8 ;
%         lim = 1e4 ;
%             while abs( ratio ) >= tol
%                 ratio = fun(chi,C,S)/fprime(chi,C,S) ;
%                 chi = chi - ratio ;
%                 ii = ii + 1 ;
%                 z = alpha*chi^2 ; % New z
%                 [ S , C ] = StumpfCalc( z , denomS , denomC ) ;
%                     if ii > lim
%                         error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
%                     end
%             end
% 
%         
% %     % Vallado
% %     bot = @(chi,C,S) chi^2*C + ( dot( r0 , v0 )/sqrt( mu ) )*chi*( 1 - z*S ) + r0mag*( 1 - z*C ) ;
% %     top = @(chi,C,S) sqrt( mu )*dt - chi^3*S - ( dot( r0 , v0 )/sqrt( mu ) )*chi^2*C - r0mag*chi*( 1 - z*S ) ; 
% %     
% %     
% %         % Newtons for UV
% %         ii = 1 ;
% %         ratio = 1 ;
% %         tol = 1e-8 ;
% %         lim = 1e4 ;
% %             while abs( ratio ) >= tol
% %                 ratio = top(chi,C,S)/bot(chi,C,S) ;
% %                 chi = chi + ratio ;
% %                 ii = ii + 1 ;
% %                 z = alpha*chi^2 ; % New z
% %                 [ S , C ] = StumpfCalc( z , denomS , denomC ) ;
% %                     if ii > lim
% %                         error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
% %                     end
% %             end
%     
% 
%    
%     
%     % Lagrange variables
%     
%     % new r
%     f = 1 - ( chi^2/r0mag )*C ;
%     g = dt - ( 1/sqrt( mu ) )*chi^3*S ;
%     r = f*r0 + g*v0 ; % new position
%     rmag = norm( r ) ; % radius
%     
%     % new v
%     fdot = ( sqrt(mu) / ( r0mag*rmag ) )*( z*chi*S - chi ) ; %(alpha*chi^3*S-chi)
%     gdot = 1 - ( chi^2/rmag )*C ; 
%     v = fdot*r0 + gdot*v0 ; % new velocity
%     
%     
%     
% 
%     
% end

function [ r1 , v1 ] = NewStateUV( r0 , v0 , dt )
mu = 398600 ;
% Find new position and velocity at some time in orbit by lagrange
% variables
    coes = state2coes( [ r0 ; v0 ] , mu ) ;
    ta = coes( 6 ) ;
    % Initial values
    r0mag = norm( r0 ) ; % radius in km
    v0mag = norm( v0 ) ; % speed in km/s
    vr0 = dot( r0 , v0 ) / r0mag ; % radial speed (km/s)   
    alpha = (2/r0mag) - ( v0mag^2 / mu ) ; % reciporical of semimajor axis
    if alpha > 0        % circle or ellipse
        h = cross( r0 , v0 ) ;
        p = norm( h )^2/mu ;
        s = 2*acot( 3*sqrt( mu / p^3 )*dt ) ;
        w = atan( nthroot( tan( s ) , 3 ) ) ;
        chi = sqrt( p )*2*cot( 2*w ) ;
    elseif alpha < 0    % hyperbola
        a = 1/alpha ;
        chi = dt*sqrt( -a )*log( ( -2*mu*alpha*dt )/( dot( r0 , v0 ) + dt*sqrt( -mu*a )*( 1 - r0mag*alpha ) ) ) ;  
        if chi > 0 && dt > 0
            chi = chi ;
        elseif chi > 0 && dt < 0
            chi = (-dt)*sqrt( -a )*log( ( -2*mu*alpha*dt )/( dot( r0 , v0 ) + (dt)*sqrt( -mu*a )*( 1 - r0mag*alpha ) ) ) ;  
        elseif chi < 0 && dt < 0 
            chi = chi ;
        elseif chi < 0 && dt > 0 
            chi = (-dt)*sqrt( -a )*log( ( -2*mu*alpha*dt )/( dot( r0 , v0 ) + (-dt)*sqrt( -mu*a )*( 1 - r0mag*alpha ) ) ) ;  
        end
    end
    
    % iterate to find chi
    var = 1 ;
    tol = 1e-8 ;
    iteration = 0 ;
    lim = 1e4 ;
    chiv = [ chi , chi ] ;
    while var > tol
        chiv(1) = chiv(2) ;
        chi = chiv(1) ;
        z = chi^2 * alpha ;
        [ c3 , c2 ] = StumpfCalc( z ) ;
        rmag = chi^2*c2 + ( dot( r0 , v0 )/sqrt( mu ) )*chi*( 1 - z*c3 ) + r0mag*( 1 - z*c2 ) ;
        dchi = ( sqrt( mu )*dt - chi^3*c3 - ( dot( r0 , v0 )/sqrt( mu ) )*chi^2*c2 - r0mag*chi*( 1 - z*c3 ) )/rmag ;
        chiv(2) = chiv(1) + dchi ;
        var = abs( chiv(2) - chiv(1) ) ;
        if iteration > lim
            error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
        end
    end
        
    f = 1 - ( chi^2/r0mag )*c2 ;
    g = dt - ( chi^3/sqrt( mu ) )*c3 ;
    gdot = 1 - ( chi^2/rmag )*c2 ;
    fdot = ( sqrt( mu )/( rmag*r0mag ) )*chi*( z*c3 - 1 )  ;
    
    r1 = f*r0 + g*v0 ;
    v1 = fdot*r0 + gdot*v0 ;
end