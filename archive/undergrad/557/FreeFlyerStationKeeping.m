clear ; close all ; clc ;
mu = 4904.8695 ;
Rmoon = 1737 ;
load( 'FreeFlyerStationKeep.mat' )
for ii = 1:length( days )-1
	deltadays(ii) = days(ii+1) - days(ii)  ;
end

I = find( deltadays >= 1 ) ;
BurnTypeUsed = BurnType(I) ;
k = 5 ;
ktgt = k ;
kint = k ;
atgt = 7500 ;
for ii = 1:length( BurnTypeUsed )
    if BurnTypeUsed(ii) == "Phase"
        dvC = CircularManeuver( atgt , SMA(I(ii)) , ECC(I(ii)) , mu ) ;
        [ ~ , ~ , dvP , ~ , ~ , ~ ] = CircularCoplanarPhasing( atgt , 10*(pi/180) , ktgt , kint , mu , Rmoon ) ;
        dv(ii) = dvP + dvC ;
    elseif BurnTypeUsed(ii) == "Circularize"
        dv(ii) = CircularManeuver( atgt , SMA(I(ii)) , ECC(I(ii)) , mu ) ;
    end
    dvTot(ii) = sum(dv(1:ii)) ;
end
figure
plot( Inclinations(I,:) )
figure
plot( RAAN(I,:) )
figure
plot( dvTot*1e3 )
figure
plot( dv*1e3 )
% figure
% plot( SMA(I,:) )
% figure
% plot( ECC(I,:) )
% figure
% plot( AOL(I,:) )

function [ tP , dtP , dvP , aP , rpP , collision ] = CircularCoplanarPhasing( atgt , phang , ktgt , kint , mu , Rcb )
    % dtP is the the time to complete phasing maneuver
    % dvP is the delta-v to go from circular to circular orbits
    %       The delta-v to go from a phasing orbit to a circular orbit
    %       is half of this
    % aP is semimajor axis of phasing orbit
    % rpP is the radius of periapsis
    % collision checks for collisions with the moon
    %       1 if a collision occurs
    %       0 if no collision
    
    ntgt = sqrt( mu / atgt^3 ) ;
    tP = ( ktgt*2*pi + phang )/ntgt ;
    aP = ( mu*( tP/( kint*2*pi ) )^2 )^(1/3) ;
    if aP < atgt
        rpP = 2*aP - atgt ;
        if rpP <= Rcb
            collision = 1 ;
        else
            collision = 0 ;
        end
    else
        rpP = atgt ;
        collision = 0 ;
    end
    dvP = 2*abs( sqrt( ( 2*mu )/( atgt ) - ( mu )/( aP ) ) - sqrt( mu / atgt ) ) ;
    dtP = ( 2*pi*ktgt ) / ( sqrt( mu / aP^3 ) ) ;
end

function dv = CircularManeuver( atgt , achase , echase , mu )
%     r = [ atgt ; 0 ] ;
    vf = [ 0 ; sqrt( mu / atgt ) ] ;
    vinit = sqrt( mu*( (2/atgt) - (1/achase) ) ) ;
    p = achase*( 1 - echase^2 ) ;
    fpa = acos( sqrt( mu*p )/( atgt*vinit ) ) ;
    vi = [ vinit*sin( fpa ) ; vinit*cos( fpa ) ] ;
    dv = norm( vf - vi ) ;
    
end