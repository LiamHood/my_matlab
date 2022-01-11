clear ; close all ; clc ;

load( 'AOmassData.mat' ) ;

sm(1,:) = data( 1:5 , 1 )' ;
sm(2,:) = data( 1:5 , 2 )' ;
sm(3,:) = data( 1:5 , 3 )' ;
sm(4,:) = data( 1:5 , 4 )' ;

em(1,:) = data( 7:11 , 1 )' ;
em(2,:) = data( 7:11 , 2 )' ;
em(3,:) = data( 7:11 , 3 )' ;
em(4,:) = data( 7:11 , 4 )' ;

for ii = 1:4
    sigs(ii) = std( sm(ii,:) ) ;
    sige(ii) = std( em(ii,:) ) ;
end

for ii = 1:4
    avgs(ii) = mean( sm(ii,:) ) ;
    avge(ii) = mean( em(ii,:) ) ;    
end

delta = -avgs+avge ;
accs = sigs./delta*100 ;
acce = sige./delta*100 ;

