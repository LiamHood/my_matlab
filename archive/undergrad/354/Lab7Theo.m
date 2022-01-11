clear ;
d = .007 ; % diameter of rod
L = .4 ; % length of rod
mb = .04128 ; % mass of block
mr = .09630 ; % mass of steel rod
E = 200e9 ; % Stiffness of steel
A = pi*(d/2)^2 ;

k = A*E/L ;
fn = (1/(2*pi))*sqrt( k/mb ) ;
