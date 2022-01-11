function [ dx ] = DoublePendulum( t , x , L , p , m , B , g )
    dx = [ x(2) ; -(B/(2*m*L^2))*x(2) + (g*p/L^2)*sin(x(1)) ] ;
end