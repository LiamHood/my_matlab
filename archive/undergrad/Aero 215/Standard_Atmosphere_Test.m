%Running function to compare to hand calculations and official numbers
%My hand calculated data and matlab data are essentially the same but the
%official number that I found diverge from my data at the higher altitudes
    [ T , P , rho ] = stdatm_HOOD_LIAM( 4500/3.28084 ); %Calculating at 4500 ft
    T1 = T;
    P1 = P;
    rho1 = rho;

    [ T , P , rho ] = stdatm_HOOD_LIAM( 32000/3.28084 );%Calculating at 32000 ft
    T2 = T;
    P2 = P;
    rho2 = rho;

    [ T , P , rho ] = stdatm_HOOD_LIAM( 42000/3.28084 );%Calculating at 42000 ft
    T3 = T;
    P3 = P;
    rho3 = rho;

    [ T , P , rho ] = stdatm_HOOD_LIAM( 82000/3.28084 );%Calculating at 82000 ft
    T4 = T;
    P4 = P;
    rho4 = rho;
    
    TM = [ T1 , T2 , T3 , T4 ]; %Matlab calculated values for temperature into a matrix
    PM = [ P1 , P2 , P3 , P4 ]; %Matlab calculated values for pressure into a matrix
    rhoM = [ rho1 , rho2 , rho3 , rho4 ]; %Matlab calculated values for density into a matrix
    column_labelM = [ "4500 ft" , "32000 ft" , "42000 ft" , "82000 ft" ]; %Altitude labels for the columns
    dataM = [ column_labelM ; TM ; PM ; rhoM ]; %All matlab data as a matrix
    labelM = [ "Altitude" ; "Matlab Temperature (K)" ; "Matlab Pressure (kPa)" ; "Matlab Density (kg/m^3)" ]; %labels for rows
    outputM = [ labelM , dataM ]; %creating output for matlab values with labels
    
    TH = [ 279.24 , 224.76 , 216.66 , 216.66 ]; %Hand calculated information
    PH = [ 85.904 , 27.467 , 17.051 , 2.495 ];
    rhoH = [ 1.0717 , .4258 , .2742 , .0401 ];
    column_labelH = [ "4500 ft" , "32000 ft" , "42000 ft" , "82000 ft" ];
    dataH = [ column_labelH ; TH ; PH ; rhoH ];
    labelH = [ "Altitude" ; "Hand Calculated Temperature (K)" ; "Hand Calculated Pressure (kPa)" ; "Hand Calculated Density (kg/m^3)" ];
    outputH = [ labelH , dataH ];
    
    TO = [ 279.24 , 224.86 , 216.66 , 216.66 ]; %Official information
    PO = [ 85.899 , 27.511 , 17.104 , 2.530 ];
    rhoO = [ 1.0717 , .4262 , .2750 , .0407 ];
    column_labelO = [ "4500 ft" , "32000 ft" , "42000 ft" , "82000 ft" ];
    dataO = [ column_labelO ; TO ; PO ; rhoO ];
    labelO = [ "Altitude" ; "Hand Calculated Temperature (K)" ; "Hand Calculated Pressure (kPa)" ; "Hand Calculated Density (kg/m^3)" ];
    outputO = [ labelO , dataO ];
    
    output = [ outputM ; outputH ; outputO ]; %Putting all data into a single matrix
    xlswrite( 'Standard_Atmosphere_Test.xlsx' , output); %Putting matrix into excel
