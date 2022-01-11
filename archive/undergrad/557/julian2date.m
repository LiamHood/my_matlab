function UTC = julian2date( JD )
    Lmonth = [ 31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 , 31 ] ;
    T1900 = ( JD - 2415019.5 )/365.25 ;
    year = 1900 + floor( T1900 ) ;
    leapyrs = floor( ( year - 1900 - 1 )*( 0.25 ) ) ;
    days = ( JD - 2415019.5 ) - (( year - 1900 )*( 365.0 )+( leapyrs )) ;
    if days < 1.0 
        year = year - 1 ;
        leapyrs = floor( ( year - 1900 - 1 )*( 0.25 ) ) ;
        days = ( JD - 2415019.5 ) - (( year - 1900 )*( 365.0 )+( leapyrs )) ;
    end
%     year_mod4 = year - floor( year/4 )*4 ;
    if mod( year , 4 ) == 0 
        Lmonth(2) = 29 ;
    end
    dayofyr = floor( days ) ;
    daysinmonths = 0 ;
    ii = 1 ;
    while sum( Lmonth(1:ii) ) < dayofyr 
        daysinmonths = sum( Lmonth(1:ii) ) ;
        ii = ii + 1 ;
    end
    month = ii ;
    day = dayofyr - daysinmonths ;
    tau = ( days - dayofyr )*24 ;
    h = floor( tau ) ;
    min = floor( ( tau - h )*60 ) ;
    s = ( tau - h - ( min/60 ) )*3600 ;
    
    UTC = [ year , month , day , h , min , s ] ;
    
end