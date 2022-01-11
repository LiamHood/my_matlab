function [ r , v ] = coes2state(  h , ecc , theta , RAAN , omega , inc , mu )
    r_peri = (h^2/mu) * ( 1/( 1 + ecc*cosd(theta) ) ) * [ cosd( theta ) ; sind( theta ) ; 0 ] ;
    v_peri = (mu/h) * [ -sind( theta ) ; ecc+cosd(theta) ; 0 ] ;

    d2r = pi/180 ;
    RAAN = d2r*RAAN ;
    omega = d2r*omega ;
    inc = d2r*inc ;
    Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
    Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
    Q(1,3) = sin(RAAN)*sin(inc) ;
    Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
    Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
    Q(2,3) = -cos(RAAN)*sin(inc) ;
    Q(3,1) = sin(inc)*sin(omega) ;
    Q(3,2) = sin(inc)*cos(omega) ;
    Q(3,3) = cos(inc) ;

    r = Q*r_peri ;
    v = Q*v_peri ;
end