function [rootH] = halley( gn , TOL , f , f1 , f2)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
x = gn ;
ii = 1 ;
for jj = 1:length(x)
    while abs( 0 - f( x(ii,jj) ) ) > TOL
        x( ii+1 , jj ) = x( ii , jj ) - (2*f(x(ii,jj))*f1(x(ii,jj)))/(2*(f1(x(ii,jj)))^2-f(x(ii,jj))*f2(x(ii,jj))) ;
        ii = ii + 1 ;
    end
    rootH(jj) = x(ii,jj) ;
end
end

