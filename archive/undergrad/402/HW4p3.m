clear ; close all ; clc ;

mdot = .05 ;
eps = 20 ;
gamma(1) = 1.4 ;
tc(1) = 300 ;
R(1) = 297 ;
gamma(2) = 1.25 ;
tc(2) = 1300 ;
R(2) = 290 ;
gamma(3) = 1.2 ;
tc(3) = 1650 ;
R(3) = 260 ;

cstar = sqrt( gamma.*R.*tc )./( gamma.*sqrt((2./(gamma+1)).^( (gamma+1)./(gamma-1) ) ) ) ;
PRlin = linspace( 2 , 400 , 1e4 ) ;
eps = ( 2/(gamma(1)+1))^(1/(gamma(1)-1)).*(PRlin).^(1/gamma(1)).*((gamma(1)+1)./(gamma(1)-1).*(1-(PRlin).^((1-gamma(1))/gamma(1)))).^(-.5) ;
index = find( eps <= 20.001 & eps >= 19.999 ) ;
PR(1) = PRlin(index) ;
eps = ( 2/( gamma(2)+1) )^( 1/( gamma(2)-1 )).*(PRlin).^(1/gamma(2)).*((gamma(2)+1)./(gamma(2)-1).*(1-(PRlin).^((1-gamma(2))/gamma(2)))).^(-.5) ;
index = find( eps <= 20.001 & eps >= 19.999 ) ;
PR(2) = PRlin(index) ;
eps = ( 2/( gamma(3)+1) )^( 1/( gamma(3)-1 )).*(PRlin).^(1/gamma(3)).*((gamma(3)+1)./(gamma(3)-1).*(1-(PRlin).^((1-gamma(3))/gamma(3)))).^(-.5) ;
index = find( eps <= 20.001 & eps >= 19.998 ) ;
 PR(3) = PRlin(index) ;
epsAct = 20 ;
cf = sqrt( ((2.*gamma.^2)./(gamma-1)).*(2./(gamma+1)).^((gamma+1)./(gamma-1)).*(1-(1./PR).^((gamma-1)./gamma)) ) + (1./PR).*epsAct ;
c = cf.*cstar ;
Isp = c ./ 9.81 ;
thrust = c.*mdot ;

properties = [ "Pressure Ratio" ; "Specific Impulse" ; "C*" ; "C_f" ; "Thrust" ] ;
ColumnNames = [ "Properties" , "Nitrogen" , "Hydrogen Peroxide" , " Hydrazine"] ;
Values = [ PR ; Isp ; cstar ; cf ; thrust ] ;
Output = [ ColumnNames ; properties , Values ]  ;
disp( Output )

% From CEA 
PRCEA = [ 284.05 , 210.81 ] ;
IspCEA = [ 1742.9/9.81 , 3200/9.81 ] ;
cfCEA = [ 1.6796 , 1.7281 ] ;
cstar = [ 1037.7 , 1851.9 ] ;

addedTitles = [ "CEA Hydrogen Peroxide" , "CEA MMH" ] ;
addedData = [ PRCEA ; IspCEA ; cstar ; cfCEA ; "???" , "???" ] ;
Output = [ Output , [ addedTitles ; addedData ] ] ;
disp( Output )
