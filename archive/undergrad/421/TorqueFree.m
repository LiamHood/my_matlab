function [ wdot ] = TorqueFree( t , w , Ix , Iy , Iz ) 
wdot = zeros( 3,1 ) ;
wdot(1) = -( Iz - Iy )*w(2)*w(3)/Ix ;
wdot(2) = -( Ix - Iz )*w(1)*w(3)/Iy ;
wdot(3) = -( Iy - Ix )*w(2)*w(1)/Iz ;
end