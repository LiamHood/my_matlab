%Liam Hood File Manipulation Lab
    
    Space_Shuttle = 'Aero_215_Lab_5_data.xlsx' ; %Assigns a variable name to the file name

    A = xlsread(Space_Shuttle , 'Space_Shuttle_Flight1' , 'A3:A14'); %Calls time for first flight and assigns to a variable
    B = xlsread(Space_Shuttle , 'Space_Shuttle_Flight1' , 'B3:B14'); %Calls altitude for first flight and assigns to a variable

    A2 = xlsread(Space_Shuttle , 'Flight2' , 'B1:M1'); %Calls time for second flight and assigns to a variable
    B2 = xlsread(Space_Shuttle , 'Flight2' , 'B2:M2'); %Calls altitude for second flight and assigns to a variable

    A3 = xlsread(Space_Shuttle , 'Flight3' , 'A1:A13'); %Calls time for third flight and assigns to a variable
    B3 = xlsread(Space_Shuttle , 'Flight3' , 'B1:B13'); %Calls altitude for third flight and assigns to a variable

    %Plotting each flight on same plot
    figure
    plot( A , B  ) %first flight plot
    hold on %plots second flight on same plot
    plot( A2 , B2 , '--' )
    hold on %plots third flight on same plot
    plot( A3 , B3 , ':k' )
    
    %Formats plot
    title( 'Altitude vs Time')
    legend( 'Flight 1' , 'Flight 2' , 'Flight 3' )
    xlabel( 'Time (s)' )
    ylabel( 'Altitude (ft)' )

    %Plotting each flight on its own plot
    figure
    subplot( 2 , 2 , 1 ) 
    plot ( A , B ) %first flight
    title( 'Altitude vs Time')
    xlabel( 'Time (s)' )
    ylabel( 'Altitude (ft)' )

    subplot( 2 , 2 , 2 )
    plot ( A2 , B2 ) %second flight
    title( 'Altitude vs Time')
    xlabel( 'Time (s)' )
    ylabel( 'Altitude (ft)' )

    subplot ( 2 , 2 , 3 )
    plot ( A3 , B3 ) %third flight
    title( 'Altitude vs Time')
    xlabel( 'Time (s)' )
    ylabel( 'Altitude (ft)' )

%Adding a new sheet with all the information
    %Changes all data into column vectors of the same length
    Ar = [0;A] ;
    Br = [0;B] ;
    A2r = [0; A2'];
    B2r = [0; B2'];
   
    
    AllData = [ Ar , Br , A2r , B2r , A3 , B3 ]; %all data becomes a single matrixj
    data = num2cell( AllData ); %converts matrix to cells
    labels = [ "Time 1 (s)  " , "Altitude 1 (ft)  " , "Time 2  " , "Altitude 2  " , "Time 3  " , "Altitude 3  " ]; %Creates labels for each column
    xlswrite( Space_Shuttle , data , 'AllFlights' , 'A2'); %Puts data into excel
    xlswrite( Space_Shuttle , labels , 'AllFlights' , 'A1'); %Puts labels on the data
    
    
