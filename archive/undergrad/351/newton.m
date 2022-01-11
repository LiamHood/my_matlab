function [ x ] = newton( x0 , f , fprime , tol , lim )
% Uses Newtons Method to find x given initial guess x0, function f,
% derivative fprime, tolerance tol, and limit on the iterations lim

    x = x0 ;
    ii = 1 ;
    ratio = 1 ;
        while abs(ratio(ii)) >= tol
            ratio(ii+1) = f(x(ii))/fprime(x(ii)) ;
            x(ii+1) = x(ii) - ratio(ii+1) ;
            ii = ii + 1 ;
                if ii > lim
                    error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
                end
        end
x = x( ii ) ;
end