%% Constellation Element Trade 35degree
% add forcemodel
clear ; close all ; clc ;

% results = xlsread( 'CTR_35d_2.csv' ) ;
% save( 'CTR_35d.mat' ) ;

load( 'CTR_35d.mat' ) ;
aR = results( : , 1 ) ;
iR = results( : , 2 ) ;
TR = results( : , 3 ) ;
PR = results( : , 4 ) ;
COVaR = results( : , 7 ) ;
COVsR = results( : , 8 ) ;
dvmaxR = results( : , 9 ) ;
dvavgR = results( : , 10 ) ;
P = [ 3 , 4 , 6 , 8 , 10 , 12 , 14 , 16 ] ;

for ii = 1:length(P)
        index = find( PR == P(ii) ) ;
        a(:,ii) = aR( index ) ;
        i(:,ii) = iR( index ) ;
        T(:,ii) = TR( index ) ;
        Ps(:,ii) = PR( index ) ;
        dvmax(:,ii) = dvmaxR( index ) ;
        dvavg(:,ii) = dvavgR( index ) ;
        a(:,ii) = aR( index ) ;
        i(:,ii) = iR( index ) ;
        T(:,ii) = TR( index ) ;
        COVa(:,ii) = COVaR( index ) ;
        COVs(:,ii) = COVsR( index ) ;
        dvmax(:,ii) = dvmaxR( index ) ;
        dvavg(:,ii) = dvavgR( index ) ;
end
aLabel = a(:,1) ;

figure
hold on
for ii = 1:length(P)
    plot( a(:,ii) , T(:,ii) )
end
hold off
xlabel( 'Semimajor Axis [km]' )
ylabel( 'Number of Satellites' )
legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' ) ;
% 
figure
hold on
for ii = 1:length(P)
    plot( a(:,ii) , dvmax(:,ii) )
end
hold off
title( 'Satellite requiring the most Delta V' )
xlabel( 'Semimajor Axis [km]' )
ylabel( 'Delta V per day [km/s/day]' )
legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' ) ;

figure
hold on
for ii = 1:length(P)
    plot( a(:,ii) , COVa(:,ii) )
end
hold off
xlabel( 'Semimajor Axis [km]' )
ylabel( 'COV' )
legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' ) ;

Tmap = mapHL( T , 9 ) ;
Pmap = mapHL( Ps , 5 ) ;
dvmap = mapHL( dvmax , 7 ) ;
COVmap = mapHL( COVa , 1 ) ;
intersection = pinter( P ) ;
D = size( a ) ;
intersectionFull = intersection ;
for ii = 2:D(1)
    intersectionFull = [ intersectionFull ; intersection ] ;
end
imap = mapHL( intersectionFull , 3 ) ;

score = Tmap + Pmap + dvmap + COVmap + imap ;
figure
surf( P , aLabel , score )
xlabel( 'Planes' )
ylabel( 'Semimajor Axis [km]' )
zlabel( 'Score in Trade Study' )





function Amap = mapHL( A , weight )
maximum = max( max( A ) ) ;
offset = min( min( A ) ) ;
D = size( A ) ;
    for ii = 1:D(1)
        for jj = 1:D(2) 
            Amap(ii,jj) = weight*( 10 - floor( 10*( A(ii,jj) - offset )/( maximum - offset ) ) ) ;
        end
    end
end

function inter = pinter( P )
ending = length(P) ;
    for jj = 1:ending
        interv = 0 ;
        for ii = 1:P(jj)-1
            interv = interv + ii ;
        end
        inter(jj) = interv*2 ;
    end
end