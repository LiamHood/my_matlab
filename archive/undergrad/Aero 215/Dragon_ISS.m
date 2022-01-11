%% Aero 215 Midterm 2
    %Liam Hood
    %Dragon capsule rendevous with the ISS
    clear all;
    clc;
    

%% Given and Initial Values
    r_circ = 10000 ; %radius of the circular orbit that both S/C begin in (km)
    r_e = 6378 ; %radius of earth (km)
    mu = 398600 ; %km^3/s^2
    ta_lead = 90 ; %degrees of true anomaly dragon is behind ISS
    
%% Delta-V for different levels of patience
  n = 1 ; %number of revolution
  while n < 10
      n = n + 1 ;
      [ delta_v ] = rendevous( r_circ , r_e , mu , ta_lead , n );
      plot( n , delta_v , '.k' )  ; hold on
      xlabel( 'Number of Revolutions' )
      ylabel( 'Delta-V in km/s' )
      title( 'Delta-V vs. Number of Revolutions' )
      
  end