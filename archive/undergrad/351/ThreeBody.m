function [ dstate ] = ThreeBody( t , state , m ) 
G = 6.67259e-20 ;
dstate = zeros( 18 , 1 ) ;
t
% 1st body state
r1 = state( 1:3 ) ;
v1 = state( 4:6 ) ;

% 2nd body state
r2 = state( 7:9 ) ;
v2 = state( 10:12 ) ;

% 3rd body state
r3 = state( 13:15 ) ;
v3 = state( 16:18 ) ;

% cube of radial differences
R12 = (norm( r2 - r1 ))^3 ;
R13 = (norm( r3 - r1 ))^3 ;
R23 = (norm( r3 - r2 ))^3 ;

% derivative of positions is velocity
dstate( 1:3 ) = v1 ;
dstate( 7:9 ) = v2 ;
dstate( 13:15 ) = v3 ;

% accelerations
dstate(4) = G*m(2)*(r2(1)-r1(1))/R12 + G*m(3)*(r3(1)-r1(1))/R13 ;
dstate(5) = G*m(2)*(r2(2)-r1(2))/R12 + G*m(3)*(r3(2)-r1(2))/R13 ;
dstate(6) = G*m(2)*(r2(3)-r1(3))/R12 + G*m(3)*(r3(3)-r1(3))/R13 ;

dstate(10) = G*m(1)*(r1(1)-r2(1))/R12 + G*m(3)*(r3(1)-r2(1))/R23 ;
dstate(11) = G*m(1)*(r1(2)-r2(2))/R12 + G*m(3)*(r3(2)-r2(2))/R23 ;
dstate(12) = G*m(1)*(r1(3)-r2(3))/R12 + G*m(3)*(r3(3)-r2(3))/R23 ;

dstate(16) = G*m(1)*(r1(1)-r3(1))/R13 + G*m(2)*(r2(1)-r3(1))/R23 ;
dstate(17) = G*m(1)*(r1(2)-r3(2))/R13 + G*m(2)*(r2(2)-r3(2))/R23 ;
dstate(18) = G*m(1)*(r1(3)-r3(3))/R13 + G*m(2)*(r2(3)-r3(3))/R23 ;

end