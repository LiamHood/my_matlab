clear ; close all ; clc ;
Thrust = 40 ; % N
Itot = 4e4 ; % N*s
g = 9.81 ;
Ptank = 20e6 ; % Mpa
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
epsa(1) = 0 ;
epsa(3) = 0 ;

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
index = 577 ;
epsa(2) = eps( 577 ) ;
eps = ( 2/(gammaN+1))^(1/(gammaN-1)).*(PR).^(1/gammaN).*((gammaN+1)/(gammaN-1).*(1-(PR).^((1-gammaN)/gammaN))).^(-.5) ;
figure 
    plot( eps , Isp )
    xlabel( 'Area Ratio' )
    ylabel( 'Specific Impulse' )
epsa(4) = eps( 577 ) ;
PR = PR( index ) ;
gamma = [ gammaH , gammaH, gammaN , gammaN ] ;
R = [ RH , RH , RN , RN ] ;
Tc = [ TcH , TcH , TcN , TcN ] ;
Cf = sqrt( ((2*gamma.^2)./(gamma-1)).*(2./(gamma+1)).^((gamma+1)./(gamma-1)).*(1-(1./PR).^((gamma-1)./gamma)) ) + (1./PR).*epsa ;
Cf(1) = sqrt( ((2*gamma(1).^2)./(gamma(1)-1)).*(2./(gamma(1)+1)).^((gamma(1)+1)./(gamma(1)-1)) ) ;
Cf(3) = sqrt( ((2*gamma(3).^2)./(gamma(3)-1)).*(2./(gamma(3)+1)).^((gamma(3)+1)./(gamma(3)-1)) ) ;
cstar = sqrt( gamma.*R.*Tc )./( gamma.*sqrt((2./(gamma+1)).^((gamma+1)./(gamma-1)) )) ;
c = Cf.*cstar ;
Ispa = c ./ g ;
mdot = Thrust./c ;
mp = Itot./(Ispa.*g) ;
tv(1:2) = mp(1:2)/1020 ;
tv(3:4) = mp(3:4)/(Ptank/(RN*TcN)) ;


properties = [ "Expansion Ratio" ; "Specific Impulse [s]" ; "Mass Flow Rate [kg/s]" ; "Propellant Mass [kg]" ; "Tank Volume [m^3]" ] ;
ColumnNames = [ "Properties" , "Hydrazine Ideal" , "Hydrazine Real" , "Cold Gas Ideal" , "Cold Gas Real" ] ;
Values = [ epsa ; Ispa ; mdot ; mp ; tv ] ;
Output = [ ColumnNames ; properties , Values ]  ;
Output(2,2) = "Inf." ;
Output(2,4) = "Inf." ;
disp( Output )

Ptank = linspace( 1e6 , 40e6 , 1e4 ) ;
tv = mp(4)./(Ptank/(RN*TcN)) ;
figure
plot( Ptank.*1e-6 , tv )
xlabel( 'Tank Pressure [MPa]' )
ylabel( 'Tank Volume [m^3]' )