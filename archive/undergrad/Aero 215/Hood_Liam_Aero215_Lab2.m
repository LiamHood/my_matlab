%Prime number finder
    
    %Finds largest prime number below the starting value of p
    p = 10^7;

    %The loop stops once a prime number is reached
    while isprime(p) == 0
        %p decreases by 1 each time through the loop
        p = p - 1;
    end
    disp([ 'The largest prime number below 10 million is ', num2str( p ) ]);

    
    
%Random number adding machine

    %sum is the sum of the random numbers added together
    sum = 0;
    
    %counter is counting how many numbers have been added
    counter = 0;
    
    %g is counting how many times a number greater than .75 has been
    %generated
    g = 0;
    
    %The while loop continues adding a new random number until sum passes
    %20
    while sum < 20
        %A new random number is generated each time through the loop
        r = rand;
            %Shows an alert, and counts, every instance of the random
            %number surpassing .75
            if r > .75
                g = g + 1 ;
                disp([ num2str(g), ' numbers greater than .75 have been generated'] )
            end
        %Counter increases by one each time through the loop    
        counter = counter + 1;
        %Sum adds the new r to the sum of the last loop until sum passes 20
        sum = r + sum;
    end 
    
    disp( [ 'It took ', num2str(counter ), ' random numbers to add past 20'] )
    

%Dynamic pressure with changing free-stream velocity

    %V is the free stream velocity, in m/s, bounded by the first and second terms. The
    %third term determines how many points the plot is calculated over
    V = linspace( 1 , 1000 , 1000 );

    %Rho is the density of air in kg/m^3
    rho = 1.225 ;

    %q is dynamic pressure defined by its relationship to rho and v
    q = .5 * rho * V.^2 ;
    
    %Creates a graph with V as the independent variable and q as dependent
    plot( V , q )
    
    %Labelling the graph
    xlabel( 'Free-Stream Velocity (m/s)' )
    ylabel( 'Dynamic Pressure (Pa)' )
    title( 'Dynamic Pressure vs Free-Stream Velocity')