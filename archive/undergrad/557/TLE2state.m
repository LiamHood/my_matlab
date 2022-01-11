function [ r , v , epoch ] = TLE2state( TLEtext )
    mu = 398600 ;
    d2r = pi/180 ;
%     fid = fopen( TLEtext , 'rb' ) ;
%     L1c = fscanf( fid , '%24c%' , 1 ) ;
%     L2c = fscanf( fid , '%71c%' , 1 ) ;
%     L3c = fscanf( fid , '%71c%' , 1 ) ;
%     fprintf( L1c ) ;
%     fprintf( L2c ) ;
%     fprintf([ L3c , '\n' ]) ;
%     fclose( fid ) ;
    fid = fopen( TLEtext , 'rb');
    L1 = fscanf( fid , '%24c%*s' , 1 );
    L2 = fscanf( fid , '%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d' , [1,9] );
    L3 = fscanf( fid , '%d%6d%f%f%f%f%f%f%f' , [1,8] );
    fclose( fid );
    
    epoch = L2(1,4)         ;        % days of year 
    Db    = L2(1,5)         ;        % Ballistic Coefficient
    inc   = L3(1,3)         ;        % Inclination [deg]
    RAAN  = L3(1,4)         ;        % Right Ascension of the Ascending Node [deg]
    ecc   = L3(1,5)/1e7     ;        % Eccentricity 
    w     = L3(1,6)         ;        % Argument of periapsis [deg]
    M     = L3(1,7)         ;        % Mean anomaly [deg]
    n     = L3(1,8)         ;        % Mean motion [Revs per day]
    
    a = (mu/(n*2*pi/(24*3600))^2)^(1/3);     % Semi-major axis [km]
    h = sqrt( mu*a*( 1 - ecc^2 ) ) ;
    % Calculate the eccentric anomaly using Mean anomaly
    err = 1e-10;            %Calculation Error
    E0 = M; t =1;
    itt = 0;
    while(t) 
           E =  M + ecc*sind(E0);
          if ( abs(E - E0) < err)
              t = 0;
          end
          E0 = E;
          itt = itt+1;
    end
    
    sin_ta = sind( E )*sqrt( 1 - ecc^2 )/( 1 - ecc*cosd( E ) ) ;
    cos_ta = ( cosd( E ) - ecc )/( 1 - ecc*cosd( E ) ) ;
    if sin_ta >= 0 && cos_ta >= 0
            ta = asind( sin_ta ) ;
        elseif sin_ta >= 0 && cos_ta <= 0
            ta = 180 - asind( sin_ta ) ;
        elseif sin_ta <= 0 && cos_ta <= 0
            ta = 180 - asind( sin_ta ) ;
        elseif sin_ta <= 0 && cos_ta >= 0
            ta = asind( sin_ta ) + 360 ;
    end
        
    COES = [ h , inc*d2r , ecc , RAAN*d2r , w*d2r , ta*d2r ] ;
    % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    [ r , v  ] = coes2state( COES , mu ) ;
end