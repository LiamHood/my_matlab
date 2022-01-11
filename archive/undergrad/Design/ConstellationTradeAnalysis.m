clear ; close all ; clc ;

% results = xlsread( 'CTResults2(2).csv' ) ;
% save( 'CTR2.mat' ) ;

load( 'CTR2.mat' ) ;
aR = results( : , 1 ) ;
iR = results( : , 2 ) ;
TR = results( : , 3 ) ;
PR = results( : , 4 ) ;
dvR = results( : , 9 ) ;
P = [ 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 ] ;

% index = ones( 30 , length( P ) ) ;
% for ii = 1:length(P) 
%     index(:,ii) = find( PR == (ii+1) ) ;
% end
% index2 = find( PR == 2 ) ;
% a = zeros( 30 , 9 ) ;
% i = zeros( 30 , 9 ) ;
% T = zeros( 30 , 9 ) ;
% dv = zeros( 30 , 9 ) ;
% a(2:end,1) = aR( index2 ) ;
% i(2:end,1) = iR( index2 ) ;
% T(2:end,1) = TR( index2 ) ;
% dv(2:end,1) = dvR( index2 ) ;

mp = 3 ;
for ii = 1:mp
        index = find( PR == (ii+1) ) ;
        a(:,ii) = aR( index ) ;
        i(:,ii) = iR( index ) ;
        T(:,ii) = TR( index ) ;
        dv(:,ii) = dvR( index ) ;
end

figure
hold on
for ii = 1:mp
    plot( a(:,ii) , T(:,ii) , '.' )
end
hold off
legend( '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' )


figure
hold on
for ii = 1:mp
    plot( a(:,ii) , dv(:,ii) , '.' )
end
hold off
legend( '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' )
