%% Pre Lab 9
% Aero 300
% Liam Hood

clear

%% Implement
x = [ 1 ; 0 ; -1 ; 0] ;
[ y , inv ] = dft( x ) ;
disp( 'DFT of x' )
disp( y )
disp( 'Inverse DFT' )
disp( inv )

%% Function

function [ y , inv ] = dft( x )
% Performs a discrete fourier transform of x where x is a column vector
    
    %defining variables for use in the function
    n = length( x ) ;
    w = exp(i*2*pi/n) ;

    %Creates the Fn matrix for transforming x
    for ii = 1:n
        for jj = 1:n
            Fn(ii,jj) = (1/sqrt(n)) * w^((ii-1)*(jj-1)) ;
            Fn_i(ii,jj) = (1/sqrt(n)) * w^((ii-1)*(jj-1)) ; %inverse transform matrix
        end
    end
    
    %Using the matrix to transform x
    y = Fn*x ;
    inv = Fn_i*x ;
end