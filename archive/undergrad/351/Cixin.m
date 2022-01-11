   clear ; close all ; clc ; 
   r1 = [ 4e4 , 0 , 0 ] ;
   v1 = [ 0 , 1.5 , 0 ] ;
    r2 = [ 0  , 0 , 0 ] ;
    v2 = [0 , 0 , 0 ] ;
    r3 = [  4e5 , 0 , 0 ] ;
    v3 = [ 0 , 1 , 0 ] ;
    r4 = [ -((4e4)-100) , 0 , 0 ] ;
    v4 = [ 0 , 1.5 , 0 ] ;
    state0 = [ r1 , v1 , r2 , v2 , r3 , v3 , r4 , v4 ]' ;
    m = [ 1e19 ; 6e24 ; 7e22 ; 1e3 ] ;
    time = 1e5 ;
    tspan = [ 0 time ] ;
    ooptions = odeset( 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 ) ;
    [ t , state2 ] = ode45( @NBody , tspan , state0 , ooptions , m ) ;
    
    figure    
    xlabel( 'X (km)' )
    ylabel( 'Y (km)' )
    zlabel( 'Z (km)' )
    speed = 2e0 ;
    for ii = 2:(length( state2 )/speed)
        drawnow
        hold on
        plot3( state2(1:speed*ii,1) , state2(1:speed*ii,2) , state2(1:speed*ii,3) , 'b')
%         plot3( state2(1:speed*ii,13) , state2(1:speed*ii,14) , state2(1:speed*ii,15) , 'k' )
        plot3( state2(1:speed*ii,19) , state2(1:speed*ii,20) , state2(1:speed*ii,21) , 'r' )
        plot3( state2(1:speed*ii,7) , state2(1:speed*ii,8) , state2(1:speed*ii,9) , 'k')        
        hold off
    end
