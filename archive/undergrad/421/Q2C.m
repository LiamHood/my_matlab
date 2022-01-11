function [ C ] = Q2C( eps , eta )
% epscross = crossmatrix( eps ) ;
% C = eye(3) + 2*epscross*epscross - 2*eta*epscross ;
% % ( 2*eta^2 - 1 )*
%     function [ across ] = crossmatrix( a )
%         across = [ 0 -a(3) a(2) ; ...
%                    a(3) 0 -a(1) ; ...
%                   -a(2) a(3) 0 ] ;
%     end
C = quat2rotm( [ eta , eps' ] ) ;
end