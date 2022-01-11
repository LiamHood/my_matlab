function [ lorenz ] = lorenzdefine( param )
%Creates lorenz equations. t is a number, pos is a row vector of [ x,y,z]
%values, param is [ sigma , rho , beta ] constants for lorenz
lorenz = @(t,pos) [ param(1) * ( pos(2) - pos(1) ) ; pos(1) * ( param(2) - pos(3) ) - pos(2) ; pos(1) * pos(2) - param(3) * pos(3) ] ;

end