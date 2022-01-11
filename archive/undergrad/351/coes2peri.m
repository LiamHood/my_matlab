function [ r , v ] = coes2peri( h , ecc , theta , mu ) 
% find perifocal position and velocity from COES. Theta in degrees
    r = (h^2/mu) * ( 1/( 1 + ecc*cosd(theta) ) ) * [ cosd( theta ) ; sind( theta ) ; 0 ] ;
    v = (mu/h) * [ -sind( theta ) ; ecc+cosd(theta) ; 0 ] ;
end