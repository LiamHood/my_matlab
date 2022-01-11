%Liam Hood
%Aero 215
%HW3 COE calculator

R = [ -6533 , 1570 , 42 ] ; %Position in km
V = [ -1.59 , -6.65 , 6.5 ] ; %Velocity in km/s^2

%Running function to calculate COEs 
[ a , e , inc , RAAN , aop , ta ] = RVtoCOE( R , V );

%displaying COEs
a = [ 'The semi-major axis is ' , num2str(a) , ' km' ] ;
disp( a )
e = [ 'The eccentricity is ' , num2str(e) ] ; 
disp( e )
inc = [ 'The inclination ' , num2str(inc) , ' degrees' ] ;
disp( inc )
RAAN = [ 'The right ascenion of ascending node is ' , num2str(RAAN) , ' degrees' ] ;
disp( RAAN )
aop = [ 'The argument of perigee is ' , num2str(aop) , ' degree' ] ;
disp( aop )
ta = [ 'The true anomaly ' , num2str(ta) , ' degrees' ] ;
disp( ta )