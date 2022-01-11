clear ; close all ; clc ;
 

tmB = readmatrix( 'HW6sB.txt' ) ;
tmT = readmatrix( 'HW6sT.txt' ) ;
tmL = readmatrix( 'HW6sL.txt' ) ;
tmR = readmatrix( 'HW6sR.txt' ) ;
% tmMv = readmatrix( 'HW6sMv.txt' ) ;
% tmMh = readmatrix( 'HW6sMh.txt' ) ;

dmh(:,1) = tmB(:,1) ;
mmh(:,1) = tmB(:,2) ;
dmh(:,2) = tmT(:,1) ;
mmh(:,2) = tmT(:,2) ;
dmv(:,1) = tmL(:,1) ;
mmv(:,1) = tmL(:,2) ;
dmv(:,2) = tmR(:,1) ;
mmv(:,2) = tmR(:,2) ;
% dmv(:,2) = tmMv(:,1) ;
% mmv(:,2) = tmMv(:,2) ;
% dmh(:,2) = tmMh(:,1) ;
% mmh(:,2) = tmMh(:,2) ;
for ii = 1:2
    mmv(:,ii) = mmv(:,ii)/(3*max( mmv(:,ii) )) ;
    mmh(:,ii) = mmh(:,ii)/(3*max( mmh(:,ii) )) ;
end
% mmv = -mmv./max( max( mmh ) ) ;
% mmh = -mmh./max( max( mmh ) ) ;
mmh(:,2) = mmh(:,2) + 9 ;
% mmh(:,3) = mmh(:,3) + 9 ;
mmv(:,2) = mmv(:,2) + 6 ;
% mmv(:,3) = mmv(:,3) + 6 ;

figure
axis square
hold on
for ii = 1:2
    plot( mmh(:,ii) , dmh(:,ii) )
    plot( dmv(:,ii) , mmv(:,ii) )
end
title( 'Shear Distribution Abaqus' )


% 
% figure
% plot( dmB , mmB )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Bottom Edge' )
% 
% figure
% plot( mmL , dmL )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Left Edge' )
% 
% figure
% plot( dmT , mmT )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Top Edge' )
% 
% figure
% plot( mmR , dmR )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Right Edge' )
% 
% figure
% plot( dmMh , mmMh )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Middle Horizontal' )
% 
% figure
% plot( mmMv , dmMv )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Middle Vertical' )
% 
