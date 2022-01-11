%% Lab 9
% Aero 300
% Liam Hood

clear
close all
%% Implement
 x = [ 1 , 0 , -1 , 0 , 0 , 9 ] ;
 [ y , inv ] = fftme( x ) ;
 disp( 'FFT of x' )
 disp( y )
 disp( 'Inverse FFT of x' )
 disp( inv )
 
 y_r = fft( x )'/sqrt(length(x));
 disp( 'Built in FFT' )
 disp( y_r )
 
 disp( 'My FFT gives the same answer as the built in function. But the ')
 disp( 'built in function doesn''t require x to be even in length. ' )
 

%% 

load( 'data.mat' )

%Defining for later use
interval = [ 0 max(t) ] ;
p = 100000;
tt = linspace( interval(1) , interval(2) , p ) ;

%Run Fast Fourier Transform
[ yfft , inv ] = fftme( X ) ;

%Find magnitude
for ii = 1:length(yfft)
    magy(ii) = sqrt( real(yfft(ii))^2 + imag(yfft(ii))^2 ) ;
end

%find phase angle
for ii = 1:length(yfft)
    ang(ii) = atan( imag(yfft(ii))/real(yfft(ii)) ) ;
end

%Plot magnitude
figure
plot( magy )
title( 'Magnitude of FFT' )
xlabel( 'Frequency' )
ylabel( 'Magnitude' )
figure
%Plot phase angle
plot( ang )
title( 'Phase angle of FFT' )
xlabel( 'Frequency' )
ylabel( 'Angle' )

%Full Power
fullpower = yfft( 2:L/2 )*2 ;

%Magnitude of full power
for ii = 1:length(fullpower)
    magf(ii) = sqrt( real(fullpower(ii))^2 + imag(fullpower(ii))^2 ) ;
end
%Phase angle of full power
for ii = 1:length(fullpower)
    angf(ii) = atan( imag(fullpower(ii))/real(fullpower(ii)) ) ;
end
%Plot of magnitude of full power
figure
plot( magf )
title( 'Magnitude of Fullpower' )
xlabel( 'Frequency' )
ylabel( 'Magnitude' )
%plot of angle of full power
figure
plot( angf )
title( 'Phase angle of Fullpower' )
xlabel( 'Frequency' )
ylabel( 'Angle' )

%Find main frequency
[ m , f ] = max( fullpower ) ;
disp( [ 'The MSD is ' , num2str(f) ] )

Pn = interp( yfft , interval , f , tt ) ; %Use formula 10.19
xinterp = dftinterp( interval , X , f , 1000 , p ) ; %Use filter formula from book on p.488

%Graph of data and Interpolations
figure
plot( t , X , '-' ,t , S , '-' , tt , Pn , '.' , tt , xinterp , '.' )
title( 'Spring Displacement vs. Time' )
xlabel( 'Time (s)' )
ylabel( 'Displacement (radians)' )
legend( 'No noise' , 'With noise' , '10.19' , 'Filter Function' )


%% Function

function [ y , inv ] = fftme( x )
% Performs a fast fourier transform of x where x is a column vector
    x(1,:) = x ;
    x(2,:) = conj( x ) ;
    for mm = 1:2
    xm = x(mm,:) ;
    %defining variables for use in the function
    n = length( xm ) ;
    w = exp((-i*2*pi)/n) ;
    mu = w^2 ;
    u = zeros(1,n/2) ;
    v = zeros(1,n/2) ;
    %Creates the small matrices for u and v
        for ii = 1:n/2
            for jj = 1:n/2
                m(ii,jj) = mu^((ii-1)*(jj-1)) ;
            end
        end
    %Seperates even and odd x components    
    for ii = 1:n/2
        xe(ii) = xm(-1+ii*2) ;
        xo(ii) = xm(ii*2) ;
    end
    %creates u and v
    for ii = 1:n/2
        u(ii) = m(ii,:) * xe';
        v(ii) = m(ii,:) * xo' ;
    end
    
    %First half of z
    for ii = 1:n/2
       z(ii) = u(ii) + w^(ii-1) * v(ii) ;
    end
    %second half of z
    for ii = 1:n/2
        z(ii+n/2) = u(ii) + w^(ii+n/2-1) * v(ii) ;
    end
 
    if mm == 1
        y = z'/sqrt(n) ; %FFT
    else
        inv = conj( z'/sqrt(n) ) ;%Inverse FFT
    end
    end   
end




function xinterp = dftinterp( interval , x , m , n , p )
% From the book in 10.3
c = interval(1);
d = interval(2);
t = c + (d-c)*(0:n-1)/n ;
tinterp = c + (d-c)*(0:p-1)/p ;
y = fft( x ) ;
yinterp = zeros(p,1) ;
yinterp( 1:m/2 ) = y( 1:m/2 ) ;
yinterp( m/2+1 ) = real( y( m/2+1 ) ) ;
if (m<n)
    yinterp( p-m/2+1) = yinterp( m/2+1 ) ;
end
yinterp( p-m/2+2:p ) = y( n-m/2+2:n ) ;
xinterp = real( ifft( yinterp ) ) * ( p/n ) ;

end

function Pn = interp( dft , interval , n , t )
% Uses formula 10.19 to interpolate. 
    a = real( dft ) ;
    b = imag( dft ) ;
    c = interval(1) ;
    d = interval(2) ;
    sum = 0 ;
    for kk = 2:n/2
        inside = (2*kk*pi*(t-c))/(d-c) ; 
        sum = sum + a(kk)*cos( inside ) - b(kk)*sin( inside ); %Find sum
    end
    Pn = (a(1)/sqrt(n)) + (2/sqrt(n))*sum + (a(n/2+1)/sqrt(n))*cos((n*pi*(t-c))/(d-c)) ;
end