[ r , v , epoch ] = TLE2state( 'SP3.txt' ) ;
COES = state2coes_display( [r;v] , 398600 ) ;