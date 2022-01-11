function [ wdot ] = GGAngVel( t , w , Ix , Iy , Iz , r ) 
Ip = [ Ix , 0 , 0 ; 0 , Iy , 0 ; 0 , 0 , Iz ] ;
T = (( 3*mu )/( norm(r)^5 ))*Rbcross*Ip*r ;
wdot = zeros( 3,1 ) ;
wdot(1) = (T(1)-( Iz - Iy )*w(2)*w(3))/Ix ;
wdot(2) = (T(2)-( Ix - Iz )*w(1)*w(3))/Iy ;
wdot(3) = (T(3)-( Iy - Ix )*w(2)*w(1))/Iz ;
end