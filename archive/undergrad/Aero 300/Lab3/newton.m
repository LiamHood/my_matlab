function [rootN] = newton( gn , TOL , f , f1)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
x = gn ;
ii = 1 ;
for jj = 1:length(x)
    while abs( 0 - f( x(ii,jj) ) ) > TOL
        x( ii+1 , jj ) = x( ii , jj ) - ( f ( x( ii , jj ) ) / f1( x( ii , jj )) ) ;
        ii = ii + 1 ;
    end
    rootN(jj) = x(ii,jj) ;
end
end

