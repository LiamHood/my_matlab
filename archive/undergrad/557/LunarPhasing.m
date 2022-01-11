clear ; close all ; clc ;
mu = 4904.8695 ;
Rcb = 1737 ;

n = 1000 ; % vary phase angle
m = 10 ; % vary revs
atgt = 7500 ;
Tnom = 2*pi*sqrt( atgt^3 / mu ) ;
ktgt = 1:1:m ;
kint = 1:1:m ;

% Set maximum angle for phasing
maxAngle = 360/8 ; % Currently set to angle between satellites on nominal orbit
phang = linspace( -maxAngle*(pi/180) , maxAngle*(pi/180) , n ) ;

dtP = zeros(n,m) ;
dvP = zeros(n,m) ;
aP = zeros(n,m) ;
rpP = zeros(n,m) ;
for ii = 1:n
    for jj = 1:m
        [ tP(ii,jj) , dtP(ii,jj) , dvP(ii,jj) , aP(ii,jj) , rpP(ii,jj) , collision ] = CircularCoplanarPhasing( atgt , phang(ii) , ktgt(jj) , kint(jj) , mu , Rcb ) ;
        % dtP is the the time to complete phasing maneuver
        % dvP is the delta-v to go from circular to circular orbits
        %       The delta-v to go from a phasing orbit to a circular orbit
        %       is half of this
        % aP is semimajor axis of phasing orbit
        % rpP is the radius of periapsis
        % collision checks for collisions with the moon
        %       1 if a collision occurs
        %       0 if no collision
    end
end
figure
surf( ktgt , phang*(180/pi) , dtP/3600 , 'EdgeColor' , 'none' )
ylabel( 'Phase Angle [degrees]' )
xlabel( 'Revolutions' )
zlabel( 'Time to Reposition [hours]' )

figure
surf( ktgt , phang*(180/pi) , rpP , 'EdgeColor' , 'none' )
ylabel( 'Phase Angle [degrees]' )
xlabel( 'Revolutions' )
zlabel( 'Radius of Perilune [km]' )

figure
surf( ktgt , phang*(180/pi) , dvP*1e3 , 'EdgeColor' , 'none' )
ylabel( 'Phase Angle [degrees]' )
xlabel( 'Revolutions' )
zlabel( 'Delta V [m/s]' )

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
