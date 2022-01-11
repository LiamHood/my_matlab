clear ; close all ; clc ;
Thrust = 40 ; % N
Itot = 4e4 ; % N*s
g = 9.81 ;
Ptank = 20 ; % Mpa
gammaH = 1.2 ;
TcH = 1650 ;
RH = 260  ;
gammaN = 1.4 ;
TcN = 300 ;
RN = 297 ;
epsa = zeros( 1 , 4 ) ;
Ispa = zeros( 1 , 4 ) ;
mdot = zeros( 1 , 4 ) ;
mp = zeros( 1 , 4 ) ;
tv = zeros( 1 , 4 ) ;
epsa(1) = "Inf." ;
epsa(2) = "Inf." ;

At = 1 ;
PR = linspace( 2 , 2e3 , 1e4 ) ;
eps = ( 2/(gammaH+1))^(1/(gammaH-1)).*(PR).^(1/gammaH).*((gammaH+1)/(gammaH-1).*(1-(PR).^((1-gammaH)/gammaH))).^(-.5) ;
Cf = sqrt( ((2*gammaH^2)/(gammaH-1))*(2/(gammaH+1))^((gammaH+1)/(gammaH-1)).*(1-(1./PR).^((gammaH-1)/gammaH)) ) + (1./PR).*eps;
cstar = sqrt( gammaH*RH*TcH )/( gammaH*sqrt((2/(gammaH+1))^((gammaH+1)/(gammaH-1)) )) ;
c = Cf.*cstar ;
Isp = c / g ;
    figure
    plot( PR , Cf ) 
    xlabel( 'Pressure Ratio' )
    ylabel( 'Coefficient of Thrust' )
    figure 
    plot( eps , Isp )
    xlabel( 'Area Ratio' )
    ylabel( 'Specific Impulse' )
    
% Isp = Thrust / ( mdot * g ) ;
% mp = Itot/(Isp*g) ;

properties = [ "Expansion Ratio" ; "Specific Impulse" ; "Mass Flow Rate" ; "Propellant Mass" ; "Tank Volume" ] ;
ColumnNames = [ "Properties" , "Hydrazine Ideal" , "Hydrazine Real" , "Cold Gas Ideal" , "Cold Gas Real" ] ;
Values = [ epsa ; Ispa ; mdot ; mp ; tv ] ;
Output = [ ColumnNames ; properties , Values ]  ;
disp( Output )