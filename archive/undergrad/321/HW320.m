r_21 = [ 6738 ; 3391 ; 1953 ] ;
z_21 = r_21/norm(r_21) ;
v_21 = [ -3.50 ; 4.39 ; 4.44 ] ;
[ r_21_cross ] = VecCross( r_21 ) ;
y_21 = (r_21_cross*v_21)/(norm(r_21_cross*v_21)) ;
[ y_21_cross ] = VecCross( y_21 ) ;
x_21 = y_21_cross*z_21 ;

x_11 = [ 1 ; 0 ; 0 ] ;
y_11 = [ 0 ; 1 ; 0 ] ;
z_11 = [ 0 ; 0 ; 1 ] ;

C_21 = [ x_21'*x_11 x_21'*y_11 x_21'*z_11 ; y_21'*x_11 y_21'*y_11 y_21'*z_11 ; z_21'*x_11 z_21'*y_11 z_21'*z_11 ];
disp( 'C_LVLH-ECI' )
disp(C_21)

alpha = acos( (trace(C_21)-1)/2 ) ;
a = [ C_21(3,2)-C_21(2,3) ; C_21(1,3)-C_21(3,1) ; C_21(2,1)-C_21(1,2) ] / 2*sin( alpha ) ;
disp( '/alpha' )
disp( alpha )
disp( 'a' )
disp( a )

theta = asin( -.62 ) ;
phi = asin( .74/cos(theta) ) ;
psi = asin( .61/cos(theta) ) ;

eta = .5*(1+trace(C_21))^.5 ;
epsilon = [ C_21(2,3)-C_21(3,2) ; C_21(3,1)-C_21(1,3) ; C_21(1,2)-C_21(2,1) ] / (4*eta) ;





function [ vec_cross ] = VecCross( vec )
    vec_cross = [ 0 , -vec(3) , vec(2) ; vec(3) , 0 , -vec(1) ; -vec(2) , vec(1) , 0 ];
end
