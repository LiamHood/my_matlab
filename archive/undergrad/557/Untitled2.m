tspan = [ 0 , 86400 ]*10 ;
mu = 398600 ;
[ r0 , v0 ] = coes2state( [ 56000 , 30*(pi/180) , 0.05 , 90*(pi/180) , 90*(pi/180) , 0 ] , mu ) ;
% r0 = [ 6000 ; 100 ; 0 ] ;
% v0 = [ 4 ; 4 ; 4 ] ;
[t, r, v] = TwoBody(tspan, r0, v0, mu, 1e-8 ) ;
[ tc , rc , vc] = CowellJ2J3( tspan , r0 , v0 , mu ) ;
[ te , re , ve ] = Encke( 10 , tspan , r0 , v0 , mu , 'gravityJ2J3' , 1 , 1 ) ;
figure 
hold on
plot3( r(1,:) , r(2,:) , r(3,:) )
plot3( re(1,:) , re(2,:) , re(3,:) )
plot3( rc(1,:) , rc(2,:) , rc(3,:) )

