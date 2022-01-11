clear ; close all ; clc ;

SMoons = [ "Janus" ; "Mimas" ; "Enceladus" ; "Tethys" ; "Dione" ; "Rhea" ; "Titan" ; "Hyperion" ; "Iapetus" ; "Phoebe" ] ;
a = [ 151.5 ; 185.5 ; 237.9 ; 294.6 ; 377.4 ; 527.1 ; 1222 ; 1481 ; 3561 ; 12960 ]*1e3 ;
P = [ 0.6947 ; 0.9424 ; 1.370 ; 1.888 ; 2.737 ; 4.518 ; 15.94 ; 21.28 ; 79.32 ; 550.6 ]*86400 ; 
Ktot = P.^2 ./ a.^3 ;
K = mean( Ktot ) ;
p = polyfit( a.^3 , P.^2 , 1 ) ;

figure
loglog( a(1:end-4).^3 , P(1:end-4).^2 , 'b*' )
hold on
loglog( a(1:end-4).^3 , (a(1:end-4).^3).*K  , 'r' )
xlabel( 'Semi-Major Axis Cubed [km^3]' )
ylabel( 'Orbital Period Squared [s^2]' )
hold off


