function OrbitPropogation( state0 , mu , time )
    tspan = [ 0 time ] ;
    options = [ 'AbsTol' , 1e-10 , 'RelTol' , 1e-10 ] ;
    [ t , state ] = ode45( @TwoBodyMotion , tspan , state0 , options , mu ) ;
    figure
    plot3( state(:,1) , state(:,2) , state(:,3) )
    xlabel( 'X (km)' )
    ylabel( 'Y (km)' )
    zlabel( 'Z (km)' )
end