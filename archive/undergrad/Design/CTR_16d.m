%% Constellation Element Trade 16 degree
% add forcemodel
clear ; close all ; clc ;

% results = readmatrix( 'CTR_16d.txt' ) ;
% save( 'CTR_16d.mat' ) ;

load( 'CTR_20d.mat' ) ;
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
for ii = 1:5
    plot( a(:,ii) , T(:,ii) )
end
% plot( a(19,2) , T(19,2) , 'ro' )
% plot( a(:,end) , T(:,end) )
hold off
% set( gca , 'XTick' , [] , 'YTick' , [] )
title( 'Number of Satellites Needed with Four Planes' )
xlabel( 'Semimajor Axis [km]' , 'fontsize' , 18)
ylabel( 'Number of Satellites' , 'fontsize' , 18)
title( '20 degree Masking Angle' , 'fontsize' , 18)
legend( '3' , '4' , '6' , '8' , '10' , '12' , '14' , '16' , 'fontsize' , 12) ;

figure
hold on
grid on
for ii = 1:4
    plot( a(:,ii) , dvmax(:,ii) )
end
plot( a(:,end) , dvmax(:,end) )
hold off
% set( gca , 'XTick' , [] , 'YTick' , [] )
title( 'Stationkeeping Delta V vs Semimajor Axis' )
xlabel( 'Semimajor Axis [km]' , 'fontsize' , 18 )
ylabel( 'Delta V per day [m/s/day]' , 'fontsize' , 18 )
legend( '3' , '4' , '6' , '8' , '16' , 'fontsize' , 12 , 'location' , 'northwest' ) ;
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

Tmap = mapHL( T , 10 ) ;
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
surf( P(1:5) , aLabel , score(:,1:5) )
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