clear ; close all ; clc ;
load( 'Group2_Debris.mat' )

piece = [ 38 ; data( 3:27 , 1 ) ] ;
Lc = data( 2:27 , 5 ) ;
A = data( 2:27 , 6 ) ;
angled = data( 2:27 , 7 ) ;
radius = data( 2:27 , 8 ) ;
color = data( 2:27 , 9 ) ;
color(1) = 0 ;
color(5) = 0 ;

angle = angled*(pi/180) ;
figure
polarplot(angle,radius,'*') 