clear ; close all ; clc ;
% 
% results = xlsread( 'EffectOfElevationAngle.xlsx' ) ;
% save( 'CTRea.mat' ) ;

load( 'CTRea.mat' ) ;
results = results( 2:end , : ) ;
aR = results( : , 1 ) ;
TR = results( : , 3 ) ;
PR = results( : , 4 ) ;
eaR= results( : , 6 ) ;
P = [ 4 , 8 , 12 , 16 ] ;

for ii = 1:length( P )
        index = find( PR == P(ii) ) ;
        a{ii} = aR( index ) ;
        T{ii} = TR( index ) ;
        ea{ii} = eaR( index ) ;
end

figure
axis( [ 0 , 50 , 0 , 200 ] )
hold on
for ii = 1:4
    plot( ea{ii} , T{ii} , '.' )
end
hold off
legend( '4' , '8' , '12' , '16' , 'Location' , 'northwest' )

