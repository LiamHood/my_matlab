function [ wdot ] = TorqueExists( t , w , Ix , Iy , Iz , T ) 
wdot = zeros( 3,1 ) ;
wdot(1) = (T(1)-( Iz - Iy )*w(2)*w(3))/Ix ;
wdot(2) = (T(2)-( Ix - Iz )*w(1)*w(3))/Iy ;
wdot(3) = (T(3)-( Iy - Ix )*w(2)*w(1))/Iz ;
end