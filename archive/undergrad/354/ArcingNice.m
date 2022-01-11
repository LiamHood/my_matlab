clear ; close all ; clc;
% data = readmatrix( 'ArcingData.xlsx' ) ;
% save( 'ArcingData.mat' )

load( 'ArcingData.mat' ) ;
dist = data( :,1 ) ;
torrin = data( :,2 ) ;
volt = data( :,3 ) ;
torr = data( :,4 ) ;
load( 'Group2_Paschen.mat' ) ;
torrin2 = [ data2( 2:7,2 ) ; data2( 9:15,2 ) ] ;
volt2 = [ data2( 2:7,4 ) ; data2( 9:15,4 ) ];
torr2 = data2( :,1 ) ;
load( 'Group3_Paschen.mat' ) ;
torrin3 = data3( 1:6,3 ) ;
volt3 = (1e3)*[ data3( 1:4,2 ) ; .529 ; data3( 6,2 ) ]  ;
torr3 = [ data3( 1:6,1 ) ; data3( 6,1 ) ];
load( 'Group4_Paschen.mat' ) ;
torrin4 = data4( 4:9,4 ) ;
volt4 = data4( 4:9,3 ) ;
torr4 = data4( :,2 ) ;

figure
hold on
loglog( torrin , volt , '*-' ) % s-s ( .224 in )
loglog( torrin2 , volt2 , '*-' ) % copper ( .273 in )
loglog( torrin3 , volt3 , '*-' ) % aluminum ( .22 in )
loglog( torrin4 , volt4 , '*-' ) % copper long (.96 in )
title( 'Paschen Curve for Stainless Steal Cathode' )
xlabel( 'Pressure-Distance (torr-in)' )
ylabel( 'Breakdown Potential (V)' )
legend( 'Stainless Steal Anode' , 'Copper Anode' , 'Aluminum Anode' , 'Copper Anode - Large Gap' , 'Location' , 'southwest' )
hold off