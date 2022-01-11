%Creates plots of temperature, pressure and density
for ii = 0:1000 %Plots 1000 points of data for each atmospheric condition
    h = ( 100000 / 1000 ) * ii; %Evenly spaces 1000 points between 0 and 100000m
    [ T , P , rho ] = stdatm_HOOD_LIAM( h ); %Runs each altitude through the function
   
    subplot( 2 , 2 , 1 ) %Creates temperature graph
    plot( T , h , '.k' ); hold on
    %labels graph
    xlabel( 'Temperature (K)' )
    ylabel( 'Altitude (m)' )
    title( 'Altitude vs Temperature' )
    
    subplot( 2 ,2 , 2)
    plot( P , h , '.b' ); hold on %creates pressure graph
    %labels graph
    xlabel( 'Pressure (kPa)' )
    ylabel( 'Altitude (m)' )
    title( 'Pressure vs Temperature' )
    
    subplot( 2 , 2 , 3 )
    %labels graph
    xlabel( 'Density (kg/m^3)' )
    ylabel( 'Altitude (m)' )
    title( 'Density vs Temperature' )
    plot( rho , h , '.r' ); hold on %creates density graph
    
  
end
