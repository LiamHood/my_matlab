clear ; close all ; clc ;
opts = odeset( 'Reltol' , 1e-8 , 'Abstol' , 1e-8 ) ;
tspan = [ 0 , 24*60*60*3 ] ;
r = [ 7.1366e03 ; 0 ; 0 ] ;
v = [ 0 ; -1.0956 ; 7.3927 ] ;
[ tvec , state ] = ode45( @TwoBody , tspan , [ r ; v ] , opts , 398600 ) ;
for ii = 1:length(tvec)
    t = tvec(ii) ;
    UTC = [ 2019 , 6 , 30 , 12 , 0 , 0 ] ;
    [ ~ , ~ , ~ , Day0 ] = Julian( UTC ) ;
    JDay = Day0 + t/(24*60*60) ;
    JYear = JDay/365.2422 ;
    [ rerels(:,ii) , ~ ] = planetary_state(3,JYear) ;
end

    theta = 23.4 ;
    Cs2eci = [ 1 0 0 ; 0 cosd(theta) sind(theta) ; 0 -sin(theta) cos(theta) ] ;
    rsun = - Cs2eci*rerels ;
    s = rsun/norm(rsun) ;
    for ii = 1:length( tvec )
        rsun(:,ii) = - Cs2eci*rerels(:,ii) ;
        s(:,ii) = rsun(:,ii)/norm(rsun(:,ii)) ;        
        r = state(1:3) ;
        directs(ii) = dot( s(:,ii) , r/norm(r) ) ;
    end
    
    figure
    plot( tvec , directs )