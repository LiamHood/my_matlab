clear ; close all ; clc ;
load( 'Lab7Vibe' ) ;
freq = VibeData( : , 1 ) ;
ref = VibeData( : , 2 ) ;
A1 = VibeData( : , 3 ) ;
A2 = VibeData( : , 5 ) ;
A3 = VibeData( : , 7 ) ;
cg = VibeData( : , 9 ) ;
P1 = VibeData( : , 4 ) ;
P2 = VibeData( : , 6 ) ;
P3 = VibeData( : , 8 ) ;

d = .007 ; % diameter of rod
L = .4 ; % length of rod
mb = .04128 ; % mass of block
mr = .09630 ; % mass of steel rod
E = 200e9 ; % Stiffness of steel
A = pi*(d/2)^2 ;
I = mb*L^2 + (1/3)*mr*L^2 ;
g = 9.81 ;
w = g*mr/L ;
P = g*mb ;
Q = 33.3 ;

k = A*E/L ;
fn = (1/(2*pi))*sqrt( k/(mb) ) ;
deltarodt = w*L^2/(24*E*I)*(3*L^2) ;
deltablockt = P*L^2/(6*E*I)*(2*L) ;
tdeltat = deltarodt + deltablockt ;
Mfn = Q*( w*L^2/2 + P*L ) ;
stresst = Mfn*(d/2)/I ;

a = A2*g ;
wa = (mr*a/L) ;
Pa = mb*a ;
deltarod = wa*L^2/(24*E*I)*(3*L^2) ;
deltablock = Pa*L^2/(6*E*I)*(2*L) ;
tdelta = deltarod + deltablock ;
Mfna = Q*( wa*L^2/2 + Pa*L ) ;
stress = Mfna*(d/2)/I ;

figure
plot( freq , tdelta )
figure
hold on 
plot( freq , A1 )
plot( freq , A2 )
plot( freq , A3 )
legend( 'Titanium' , 'Steel' , 'Aluminum' )
xlabel( 'Frequency (Hz)' )
ylabel( 'Acceleration (g)' )
title( 'Vibe Table' )
hold off
