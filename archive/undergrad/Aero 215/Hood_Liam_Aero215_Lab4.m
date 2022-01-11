%Cross product calculator
    %Input vectors
    H = [ 1 2 3 ] ;
    N = [ 7 8 9 ] ;

    %Using cross product function
   [ crossed ] = cross_product( H,N ) ; 
   
    %Returning the answer
    disp( 'My answer' ) ;
    disp(crossed) ;
    
    %Checking 
    actual_cross = cross( H , N ) ;
    
    disp('Matlab answer') ;
    disp(actual_cross) ;
    
   
%Velocity power plot
    %Velocity vector
    V = [ 1 2 3 4 5 ] ;
    
    %Time vector
    t = [ .2 .4 .6 .9 1.1 ] ;
    
    %Plotting the different powers of velocity against time
    plot( t , V , t , V.^2 , t , V.^3 );
    
    %Making the graph make sense
    legend( 'V vs. t' , 'V^2 vs. t' , 'V^3 vs. t' );
    title( 'Powers of Velocity vs. Time' );
    xlabel( 'time (s)' )
    ylabel( 'Velocity (ft/sec)' )
    
    
    
   
