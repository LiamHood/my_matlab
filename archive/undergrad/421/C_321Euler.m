function [ C ] = C_321Euler( roll , pitch , yaw )
% Find rotation matrix for 3-2-1 Euler rotation from yaw pitch and roll
C = zeros( 3 , 3 ) ;
C(1,:) = [ cos( pitch )*cos( yaw ) , cos( pitch )*sin( yaw ) , -sin(pitch) ] ;
C(2,1) = sin( roll )*sin( pitch )*cos( yaw ) - cos( roll )*sin( yaw ) ;
C(2,2) = sin( roll )*sin( pitch )*sin( yaw ) + cos( roll )*cos( yaw ) ;
C(2,3) = sin( roll )*cos( pitch ) ;
C(3,1) = cos( roll )*sin( pitch )*cos( yaw ) + sin( roll )*sin( yaw ) ;
C(3,2) = cos( roll )*sin( pitch )*sin( yaw ) - sin( roll )*cos( yaw ) ;
C(3,3) = cos( roll )*cos( pitch ) ;

end