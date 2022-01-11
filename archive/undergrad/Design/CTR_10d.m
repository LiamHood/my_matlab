%% Constellation Element Trade 10 degree
% add forcemodel
clear ; close all ; clc ;

% results = xlsread( 'CTR_20d_2.csv' ) ;
% save( 'CTR_0d.mat' ) ;

load( 'CTR_10d.mat' ) ;
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
grid on
for ii = 1:5 %length(P)
    plot( a(:,ii) , T(:,ii) )
end
hold off
% set( gca , 'XTick' , [] , 'YTick' , [] )
xlabel( 'Semimajor Axis [km]' , 'fontsize' , 18)
ylabel( 'Number of Satellites' , 'fontsize' , 18)
title( '10 degree Masking Angle' , 'fontsize' , 18)
legend( '3 Planes' , '4 Planes' , '6 Planes' , '8 Planes' , '10 Planes' , '12 Planes' , '14 Planes' , '16 Planes' , 'fontsize' , 12) ;

figure
hold on
grid on
for ii = 1:length(P)
    plot( a(:,ii) , dvmax(:,ii) )
end
hold off
% set( gca , 'XTick' , [] , 'YTick' , [] )
xlabel( 'Semimajor Axis [km]' , 'fontsize' , 18 )
ylabel( 'Drift of Planes' , 'fontsize' , 18 )
legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' , 'fontsize' , 12 , 'location' , 'northwest' ) ;
% 
% figure
% hold on
% for ii = 1:length(P)
%     plot( a(:,ii) , COVa(:,ii) )
% end
% hold off
% xlabel( 'Semimajor Axis [km]' )
% ylabel( 'COV' )
% legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' ) ;

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