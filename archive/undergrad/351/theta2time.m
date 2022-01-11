function [ time ] = theta2time( theta , T , ecc )
theta = theta*(pi/180) ;
n = 2*pi/T ;
E = 2*atan( sqrt((1-ecc)/(1+ecc)) * tan( theta/2 ) ) ;
time = ( E - ecc*sin(E) ) / n ;
    if time < 0
        time = T+time ;
    end
end