[t, data] = serial_reader(100) ;
save( 'Data_Distance_Speed' )
figure
plot( data(:,1) , data(:,2) )
title( 'Speed vs distance' )
xlabel( 'Distance (cm)' )
ylabel( 'Speed' )

figure
plot( t , data(:,1) , t , data(:,2) )
title( 'Speed and distance vs time' )
ylabel( 'Distance (cm) and Speed' )
xlabel( 'Time' )