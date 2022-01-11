%Lab 5
%Aero 300
%Liam Hood

clear 
close all

%% Inputs
airfoil_1 = textread( 'airfoil_1.txt' ) ;
airfoil_2 = textread( 'airfoil_2.txt' ) ;
cd = textread( 'cd_data.txt' ) ;
cl = textread( 'cl_data.txt' ) ;

%% 1
%Interpolate root of wing values
RootWingC = polyfit( airfoil_1( : , 1 ) , airfoil_1( : , 2 ) , 5 ) ;
RootWing = @(x) RootWingC(1)*x.^5 + RootWingC(2)*x.^4 + RootWingC(3)*x.^3 + RootWingC(4)*x.^2 + RootWingC(5)*x + RootWingC(6) ;

%Interpolate tip of wing values
TipWingC = polyfit( airfoil_2( : , 1 ) , airfoil_2( : , 2 ) , 5 ) ;
TipWing = @(x) TipWingC(1)*x.^5 + TipWingC(2)*x.^4 + TipWingC(3)*x.^3 + TipWingC(4)*x.^2 + TipWingC(5)*x + TipWingC(6) ;

n = 100 ;
intervalx = linspace( 0 , 1 , n ) ;
intervaly = linspace( 0 , 1 , n ) ;
for ii = 1:n
    WingHeight(:,ii) = linspace( RootWing( intervaly(ii) ) , TipWing( intervalx(ii) ) , n ) ;
end
for jj = 1:n
    for kk = 1:n
        if WingHeight( jj , kk ) < 0 
            WingHeight( jj , kk ) = 0;
        end
    end
end
% Wing( : , 1 ) = RootWing( intervaly ) ;
% Wing( : , n ) = TipWing( intervaly ) ;


surf( intervalx , intervaly , WingHeight , 'lineStyle' , 'none' ) ; hold on
surf( intervalx , intervaly , -1.*WingHeight , 'lineStyle' , 'none' ); hold on
axis( [ 0 , 1 , 0 , 1 , -.5 , .5 ] )

%% 2
f = find( cl(:,2) < .834 ) ; %Finding nice lift data that follows thin foil theory
cl_c = polyfit( cl(f,1) , cl(f,2) , 2 ) ; %Polynomial coeffecients fitting lift data best
cl_f = @(x) cl_c(1)*x.^2 + cl_c(2)*x +cl_c(3) ; %function that best fits lift data
disp([ 'The zero angle of attack lift coeffecient is ' , num2str( cl_c( 3 ) ) ]) ;
cl_ang = cl(f,1) ;%angles of attack with well behaved cl values
interval_cl = linspace( min( cl_ang ) , max( cl_ang ) , 100 ) ; 

cd_c = polyfit( cd(:,1) , cd(:,2) , 2 ) ; %Polynomial coeffecients fitting drag data best
cd_f = @(x) cd_c(1)*x.^2 + cd_c(2)*x +cd_c(3) ; %function that best fits drag data
interval_cd = linspace( min( cd(:,1) ) , max( cd(:,1) ) , 100 ) ;


figure
hold on
plot( cd(:,1) , cd(:,2) , 'r.' ) %Plot of drag data
plot( cl(:,1) , cl(:,2) , 'b.' ) %Plot of lift data
plot( interval_cl , cl_f( interval_cl ) , 'b-' ) %Plot of lift interpolation
plot( interval_cd , cd_f( interval_cd ) , 'r-' ) %Plot of lift interpolation

%Labeling graph
title( 'Coeffecients of Lift and Drag vs Angle of Attack' )
xlabel( 'Angle of attack' )
ylabel( 'Coeffecient Value' )
legend( 'Coeffecient of Drag Data' , 'Coeffecient of Lift Data' , 'Coeffecient of Drag Interpolation' , 'Coeffecient of Lift Interpolation' , 'Location' , 'northwest' ) ;


