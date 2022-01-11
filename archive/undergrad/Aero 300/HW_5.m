%% 5.1.C1
for ii = 1:12
    h(ii) = 10^(-ii) ; %Step sizes used for the following calculation
end

f = @(x) sin(x) - cos(x) ; %Original function
fprime = @(x) cos(x) + sin(x) ; %first derivative
fdp = @(x) -sin(x) + cos(x) ; %second derivative
ftp = @(x) -cos(x) - sin(x) ; %third derivative

fpEstimate = ( f(h) - f(-h) ) ./ ( 2*h ) ;

fpActual = linspace( fprime(0) , fprime(0) , 12 ) ;

eTheoreticalMax = ( ( h.^2 ) / 6 ) .* ftp(h) ; %Error if c is h
eTheoreticalMin = ( ( h.^2 ) / 6 ) .* ftp(-h) ; % Error if c is -h

eActual = fpEstimate - fpActual ; %True error

%column labels for table
ErrorTableTitles = [ "h" , "Actual Derivative Values" , "Estimated Derivative Values" , "Theoretical Max Error" , "Theoretical Min Error" , "Actual Error" ];
%data columns for table
ErrorTableData = [ h' , fpActual' , fpEstimate' , eTheoreticalMax' , eTheoreticalMin' , eActual' ];
ErrorTable = [ ErrorTableTitles ; ErrorTableData ] ; %combining titles and data into one table
disp( ErrorTable ) %display table

plot( h , eTheoreticalMax , '-r' , h , eTheoreticalMin , '-b' , h , eActual , '-k' )
title( 'Error vs Step Size' )
xlabel( 'Step Size' )
ylabel( 'Error' )
legend( 'Max Error' , 'Min Error' , 'True Error' ) 

disp( 'Actual error is withing bounds of theoretical error' )

%% 5.2.1a
%Inputs
a = 0 ; %lower bound of interval of integration
b = 4 ; %upper bound of interval of integration
f = @(x) x/( sqrt( x^2 + 9 ) ) ; %function
AntiDerivative_f = @(x) sqrt( x^2 + 9 ) ; %actual antiderivative
Integral_f = AntiDerivative_f(4) - AntiDerivative_f(0) ; %Actual integral
%Run for 16 panels
Int_m16 = ctm( f , a , b , 16 ); 
Error_16 = Integral_f - Int_m16 ;
disp([ 'Using the trapezoid rule with 16 panels gives us a value of ' , num2str( Int_m16 ) , ' with an error of ' , num2str( Error_16 ) ]) ;

%Run for 32 panels
Int_m32 = ctm( f , a , b , 32 );
Error_32 = Integral_f - Int_m32 ;
disp([ 'Using the trapezoid rule with 32 panels gives us a value of ' , num2str( Int_m32 ) , ' with an error of ' , num2str( Error_32 ) ]) ;


function [ ctm_f ] = ctm( f , a , b , m )
%Uses composite trapezoid method to calculate the integral of f over the
%interval of a to b with m panels
    h = ( b - a ) / m ; %Step size
    sum_f = 0 ; %starting value for sum 
    for ii = 1:(m-1)
        sum_f = sum_f + f( a + h*ii );  %sum in trapezoid method
    end
    ctm_f = ( h/2 ) * ( f(a) + f(b) + 2*sum_f ) ; %composite trapezoid formula
end


