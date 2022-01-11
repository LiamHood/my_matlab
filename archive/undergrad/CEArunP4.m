Ae_At = [ 1.0000 ; 1.0000 ; 2.0000 ; 5.0000 ; 10.000 ; 100.00 ] ; 
Isp = [ 1523.3 ; 1523.3 ; 2805.6 ; 3451.6 ; 3786.1 ; 4448.7 ] ;
Cf = [ 0.6569   0.6569   1.2100   1.4885   1.6328   1.9186 ] ;

figure
loglog( Ae_At , Isp ) 
xlabel( 'Exit area to throat ratio' )
ylabel( 'Isp (M/s)' )
title( 'Isp vs. Area Ratio' )


figure
loglog( Ae_At , Cf ) 
xlabel( 'Exit area to throat ratio' )
ylabel( 'Coefficient of thrust' )
title( 'C_f vs. Area Ratio' )

