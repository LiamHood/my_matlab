%Liam Hood
%Aero 300
%HW 1

c = ones( 1 , 51 ) ; %Coeffecients of p
d = 50 ; %Degree of p 
x = 1.00001 ; %The point to evaluate at
p =  nest( d , c , x ) ; %Evaluates using Horner's Methon
disp([ 'P is ' , num2str(p) ])
q = ( x^51 - 1 ) / ( x - 1 ) ; %Evaluates an equivalent equation
disp([ 'Q is ' , num2str(q) ])
error = abs( q - p ) ; %Calculates error
disp([ 'The error of p is ' , num2str(error) ])

