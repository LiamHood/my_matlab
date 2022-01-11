% PreLab 5
% Aero 300
% Liam Hood

%Inputs
x0 = [ -5 -2 4 5 ] ; %x input data
y0 = [ 4 -1 2 -5 ] ; %y input data
interval =  -5 : .1 : 5  ; %interval of evaluation
n = 3 ; %degree of polynomial to be evaluated

%Find best fit polynomial
pc = polyfit( x0 , y0 , n )  ; %Find coeffecients
p = @(x) pc(1)*(x.^3) + pc(2)*(x.^2) + pc(3)*x + pc(4) ; %Creating polynomial function

%Plotting
hold on
plot( interval , p(interval) ) %Plots interpolation function
plot( x0 , y0 , 'r.' ) %Plots original data points as red dots
title( 'Interpolation of 4 points' ) %Titles plot
xlabel( 'x' ) %labels x axis
ylabel( 'y' ) %labels y axis


disp( 'Polyfit uses QR factorization. This is order n^3 while Newtons divided difference is order n^2' )