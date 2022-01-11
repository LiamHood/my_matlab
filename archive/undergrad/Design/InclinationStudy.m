clear ; close all ; clc ;

% results = xlsread( 'CTR_IncStudy.csv' ) ;
% save( 'CTR_is.mat' ) ;

load( 'CTR_is.mat' ) ;
aR = results( : , 1 ) ;
iR = results( : , 2 ) ;
TR = results( : , 3 ) ;
PR = results( : , 4 ) ;
eaR= results( : , 6 ) ;
COVR = results( :,7 ) ;
dvR = results( : , 10 ) ;
P = [ 3 , 4 , 6 , 8 ] ;

for ii = 1:length( P )
        index = find( PR == P(ii) ) ;
        i(:,ii) = iR( index ) ;
        T(:,ii) = TR( index ) ;
        dv(:,ii) = dvR( index ) ;
end

figure
axis( [ 35 , 90 , 0 , 100 ] )
hold on
for ii = 2:4
    plot( i(:,ii) , T(:,ii) )
end
hold off
legend( '4 Planes' , '6 Planes' , '8 Planes' , 'Location' , 'northwest' )
xlabel( 'Inclination [degrees]' )
ylabel( 'Number of Satellites' )


figure
hold on
for ii = 2:4
    plot( i(:,ii) , dv(:,ii) )
end
hold off
legend( '4' , '6' , '8'  , 'Location' , 'northwest' )
xlabel( 'Inclination [degrees]' )
ylabel( 'Delta V to Stationkeeping' )

