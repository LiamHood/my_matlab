%% Given values
        s = 24 ; %wing area in m^2
        q = .27 * 10^3 ; %dynamic pressure in Pa
        CDo = .075 ; %zero-lift drag coefficient
        sa = -5 ; %minimum angle of attack in degrees
        la = 15 ; %maximum angle of attack in degrees
        l_sa = -177.5 ; %lift at min angle of attack in Newtons
        l_la = 5452 ; %lift at max angle of attack in Newtons
        
%% Part A
disp( 'Part "a" calculates minimum and maximum coeffecients of lift' )
        CL_min = l_sa / ( q * s ) ; %minimum coeffecient of lift
        CL_max = l_la / ( q * s ) ; %maximum coeffecient of lift
        
%% Part B and C
disp( 'Part "b" and "c" create a graph' )
disp( '     See Figure 1' )
    CL = linspace( CL_min , CL_max , 100 ) ; %Coefficients of lift
    k = .4; % k value
    while k <= 1.2
        
        CD = CDo + k .* CL .^ 2 ; %Coefficients of drag
        
        %Plots drag polars with dots at L/D max values
        hold on
        plot( CD , CL ) %Plots drag polars
        LD =  CL ./ CD ; %Calculates L/D for all CL and CD values
        [ LD_max , I ] = max( LD ) ; %finds maximum LD value for each curve
        plot( CD( I ) , CL ( I ) , '.k' ) %plots all the L/D max values onto their corresponding curves
        
        %Records all LD_max values in a vector for use outside the loop
        if k == .4
            LD_maxv(1) = LD_max;
        elseif k == .6
            LD_maxv(2) = LD_max;
        elseif k == .8
            LD_maxv(3) = LD_max;
        elseif k == 1.0
            LD_maxv(4) = LD_max;
        elseif k == 1.2
            LD_maxv(5) = LD_max;
        end
        k = k + .2 ;
    end
    
    %labels for the graph
    title ( 'Drag Polars' )
    xlabel( 'Coefficient of Drag' )
    ylabel( 'Coefficient of Lift' )
    legend( 'k = .4' , 'L/D max ' , 'k = .6' , 'L/D max ' , 'k = .8' , 'L/D max ' , 'k = 1.0' , 'L/D max ' , 'k = 1.2' , 'L/D max ' ) 
    
%% Part D
disp( 'Part d' )
    [ LD_max , k ] = max( LD_maxv ); %finds best of the L/D max values and associated curve
    k = .2 + .2 * k ; %converts the index of the curve that the LD_max value came from to its corresponding k value
    disp( [ 'The best L/D max is ' , num2str( LD_max ) , ' with a k value of ' , num2str( k ) ] )
    
    