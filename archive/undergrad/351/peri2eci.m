function [ r , v ] = peri2eci( r_peri , v_peri , RAAN , omega , inc )
% finds state vectors from COES and perifocal coordinates. Angles in
% degrees
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