clear ; clc ; close all ;
Pstd = 14.696 ; % Ambient air pressure
Tstd = 519.7 ; % Ambient air temperature
f2r = 458.67 ;
T1 = 67 + f2r ;
P01 = 14.7346 ;
R = 1717 ;
Cp = .24*(32.174) ;
k = 1.4 ;
rho = Pstd/( R*Tstd ) ;

A = 5.5 ; % square inches
%  data = xlsread( 'Lab2.xlsx' ) ;
 load( "data.mat" )
%  data = data( 1:end , : ) ;
 Throttled = data( 5:10 , : ) ;
 Nozzled = data( 11:13 , : ) ;
 throttle = Throttled( : , 1 ) ;
 P1 = Throttled( : , 4 ) + P01 ;
 T03 = Throttled( : , 8 ) + f2r ;
 P03 = Throttled( : , 10 ) + P01 ;
 T04 = Throttled( : , 14 ) + f2r;
 T05 = Throttled( : , 16 ) + f2r ;
 
 thrust_mv = data( 11:13 , 20:21 ) ;
 for ii = 1:3
     thrust_mvavg(ii) = mean( thrust_mv( ii , 1:2 ) ) ;
 end
 thrust = thrust_mvavg/.12 ;
 PercentOpen = [ 0 , 50 , 100 ] ;
 
 CompressorRatio = P03 ./ Pstd ;
 v = sqrt( ( P01 - P1 )*2/rho ) ;
 mdot = rho*A*v ;
 Corrected = mdot*( sqrt(T1/Tstd) / (P01/Pstd) ) ; % assuming sea level conditions
 
 PRcomp = P03/P1 ;
 h1 = Cp*T1 ;
 h2s = Cp*(( P03./P01 ).^((k-1)/k))*T1 ;
 h2 = Cp*T03 ;
 compeff = ( h2s - h1 )/( h2 - h1 ) ;
 
 
 figure
 plot( Corrected , PRcomp ) 
 xlabel( "Corrected Mass Flow [slug/s]" )
 ylabel( "Pressure Ratio" )
 
 figure 
  plot( Corrected , compeff ) 
 xlabel( "Corrected Mass Flow [slug/s]" )
 ylabel( "Compressor Efficiency" )

 figure
 plot( PercentOpen , thrust )
 xlabel( 'Openness of Nozzle' )
 ylabel( 'Thrust [lbf]' )
 
 
 
 %% Functions
 