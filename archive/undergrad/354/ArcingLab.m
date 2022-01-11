clear ; close all ; clc;
% data = readmatrix( 'ArcingData.xlsx' ) ;
% save( 'ArcingData.mat' )

load( 'ArcingData.mat' ) ;
dist = data( :,1 ) ;
torrin = data( :,2 ) ;
volt = data( :,3 ) ;
torr = data( :,4 ) ;
load( 'Group2_Paschen.mat' ) ;
torrin2 = data2( :,2 ) ;
volt2 = data2( :,4 ) ;
torr2 = data2( :,1 ) ;
load( 'Group3_Paschen.mat' ) ;
torrin3 = [ data3( 1:4,3 ) ; data3( 6,3 ) ] ;
volt3 = (1e3)*[ data3( 1:4,2 ) ; data3( 6,2 ) ]  ;
torr3 = [ data3( 1:4,1 ) ; data3( 6,1 ) ];
load( 'Group4_Paschen.mat' ) ;
torrin4 = data4( :,4 ) ;
volt4 = data4( :,3 ) ;
torr4 = data4( :,2 ) ;

figure
hold on
loglog( torrin , volt , '*' ) % s-s ( .224 in )
loglog( torrin2 , volt2 , '*' ) % copper ( .273 in )
loglog( torrin3 , volt3 , '*' ) % aluminum ( .22 in )
loglog( torrin4 , volt4 , '*' ) % copper long (.96 in )
title( 'Paschen Curve' )
xlabel( 'Pressure-Distance (torr-in)' )
ylabel( 'Breakdown Potential (V)' )
legend( 'Group1' , 'Group2' , 'Group3' , 'Group4' )

n = 3 ;
torrinches = linspace( 0 , 1 , 1e2 ) ;
P1 = polyfit( torrin , volt , n ) ;
volts1 = P1(1)*torrinches.^n + P1(2)*torrinches.^(n-1) + P1(3)*torrinches + P1(4) ;

P2 = polyfit( torrin2 , volt2 , n ) ;
volts2 = P2(1)*torrinches.^n + P2(2)*torrinches.^(n-1) + P2(3)*torrinches + P2(4) ;

P3 = polyfit( torrin3 , volt3 , n ) ;
volts3 = P3(1)*torrinches.^n + P3(2)*torrinches.^(n-1) + P3(3)*torrinches + P3(4) ;

P4 = polyfit( torrin4 , volt4 , n ) ;
volts4 = P4(1)*torrinches.^n + P4(2)*torrinches.^(n-1) + P4(3)*torrinches + P4(4) ;

loglog( torrinches , volts1 ) % s-s ( .224 in )
loglog( torrinches , volts2 ) % copper ( .273 in )
loglog( torrinches , volts3 ) % aluminum ( .22 in )
loglog( torrinches , volts4 ) % copper long (.96 in )

legend( 'Group1' , 'Group2' , 'Group3' , 'Group4' , 'Group1' , 'Group2' , 'Group3' , 'Group4' )
hold off

