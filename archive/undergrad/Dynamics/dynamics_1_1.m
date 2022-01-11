%Givens
p = -50;    %Newtons
m = 50;     %kg
g = 9.8;    %m/s^2
mus = .20;
muk = .15;


while p < 250
    p = p + .1;
    if p - m * g * sind( 15 ) < m * g * cosd( 15 ) * mus
        a = 0;
    elseif p - m * g * sind( 15 ) == m * g * cosd( 15 ) * mus
        disp([ 'Acceleration is zero, because of friction if p is ' , p , ' newtons or less' ])
    else
        a = ( p - m * g * cosd( 15 ) * muk - m * g * sind( 15 ) ) / m;
        hold on
        plot( p , a , '.b' ) 
       
    end
end 

title( 'Box acceleration vs P' )
xlabel( 'P (N)' )
ylabel( 'Box acceleration (m/s^2)' )