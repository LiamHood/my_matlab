close all
clear
T = [ 290 ; 440 ; 860 ; 670 ] ;
p_i_g = [ .11 ; 31 ; 31 ; 2.15 ] ;
R = [ 287 ; 287 ; 300 ; 300 ] ;
psia = 29.79 ;
p_i = p_i_g + psia ;
p = p_i * 6894.76 ;
%p(4) = 4.0097e+05 ;
v = ( R .* T ) ./ p ;
s = [ 1.66802 ; 2.0887 ; 2.79783 ; 2.654565 ] ;
s_calc(1) = 1.66802 ;
for ii = 2:4
    s_calc(ii) = s_calc(ii-1) + .718*log( T(ii) / T(ii-1) ) + .287*log( v(ii) / v(ii-1) ) ;
end

figure
%plot( s , T , '.-b' , s_calc , T , '.-r' , [s(1),s(5)] , [T(1),T(5)] , '--b' , [s_calc(1),s_calc(5)] , [T(1),T(5)] , '--r' )
plot( s_calc , T , '.-b' , [s_calc(1),s_calc(4)] , [T(1),T(4)] , '--b' , s , T )
axis([ 1.5 , 3 , 0 , 1000 ])
%legend( 's values from table A-22' , 's values by calculation' , 'Location' , 'southeast' )
title( 'T-s' )
xlabel( 'Specific Entropy ( kJ/(kg*K) )' )
ylabel( 'Temperature (K)' )
figure
%plot( [v(1),v(2),v(3),v(5)] , [p(1),p(2),p(3),p(5)]*(10^-3) , '.-b' , [v(5),v(1)] , [p(5),p(1)]*10^-3 , '--b' )
plot( v , p*(10^-3) , '.-b' , [v(4),v(1)] , [p(4),p(1)]*10^-3 , '--b' )
axis([ 0 , 1.5 , 0 , 450 ])
title( 'p-v' )
xlabel( 'Specific Volume (m^3/kg)' ) 
ylabel( 'Pressure (kPa)' )


u = [ 206.91 , 315.30 , 641.40 , 533.645 ] ;
h = u' + p.*v ;
h2 = [ 290.16 , 441.61 , 888.27 , 681.14 ] ;
powerv = .037969 * ((h(3)-h(4))-(h(2)-h(1)))
powerv2 = .037969 * ((h2(3)-h2(4))-(h2(2)-h2(1))) * 1000