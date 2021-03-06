function RocketAndWheelPlots( tspan , state , opts , dwrel , Td , Iw , Ir , name )

    [t,state]=ode45( @RocketAndWheel , tspan , state , opts , dwrel , Td , Iw , Ir ) ;

    figure( 'Name' , name , 'NumberTitle' , 'off' , 'Position' , [ 100 50 1000 650 ] ) ;

    subplot(2,3,1)
    plot( t , state(:,1) )
    title( 'Angular Displacemet X vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Displacement X (rad/s)' )
    
    subplot(2,3,2)
    plot( t , state(:,2) )
    title( 'Angular Displacemet Y vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Displacement Y (rad/s)' )
    
    subplot(2,3,3)
    plot( t , state(:,3) )
    title( 'Angular Displacemet Z vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Displacement Z (rad/s)' )
    
    subplot(2,3,4)
    plot(t,state(:,4))
    title( 'Angular Velocity X vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Velocity X (rad/s)' )
    
    subplot(2,3,5)
    plot(t,state(:,5))
    title( 'Angular Velocity Y vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Velocity Y (rad/s)' )
    
    subplot(2,3,6)
    plot(t,state(:,6))
    title( 'Angular Velocity Z vs Time' )
    xlabel( 'Time(s)' )
    ylabel( 'Angular Velocity Z (rad/s)' )

end


