clear

fun = @(t,y) ( y - t - 1 )^2 + 2 ;
tspan = [ 0 pi/3 ] ;
y0 = 1 ;
h = .01 ;
rTol = 10^-6 ;

[ tODE , yODE ] = ode45( fun , tspan , y0 );

[t,y] = dp45(fun, tspan, y0, h, rTol) ;
plot( t , y , '-r' , tODE , yODE , '-b' ) 

function [t,y] = dp45(fun, tspan, y0, h, rTol)
% Runge-Katta 4 with 5th order check
w(:,1) = y0 ;
z(:,1) = y0 ;
t(1) = tspan(1) ;
ii = 1 ;
jj = 1 ;
hIN = h ;
c31 = .3 ;
c32 = 3/40 ;
c33 = 9/40 ;
c41 = 44/45 ;
c42 = 56/15 ;
c43 = 32/9 ;
c51 = 8/9 ;
c52 = 19372/6561 ;
c53 = 25360/2187 ;
c54 = 64448/6561 ;
c55 = 212/729 ;
c61 = 9017/3168 ;
c62 = 355/33 ;
c63 = 46732/5247 ;
c64 = 49/176 ;
c65 = 5103/18656 ;
cz1 = 35/384 ;
cz2 = 500/1113 ;
cz3 = 125/192 ;
cz4 = 2187/6784 ;
cz5 = 11/84 ;
ce1 = 71/57600 ;
ce2 = 71/16695 ;
ce3 = 71/1920 ;
ce4 = 17253/339200 ;
ce5 = 22/525 ;
ce6 = 1/40 ;

    while t(ii) < tspan(2)
        s1(:,ii) = fun( t(ii) , z(:,ii) ) ;
        s2(:,ii) = fun( t(ii) + .2*h , z(:,ii) + .2*h*s1(:,ii) ) ;
        s3(:,ii) = fun( t(ii) + c31*h , z(:,ii) + c32*h*s1(:,ii) + c33*s2(:,ii) ) ;
        s4(:,ii) = fun( t(ii) + .8*h , z(:,ii) + c41*h*s1(:,ii) - c42*h*s2(:,ii) + c43*h*s3(:,ii) ) ;
        s5(:,ii) = fun( t(ii) + c51*h , z(:,ii) + h*( c52*s1(:,ii) - c53*s2(:,ii) + c54*s3(:,ii) - c55*s4(:,ii) ) ) ;
        s6(:,ii) = fun( t(ii) + h , z(:,ii) + h*( c61*s1(:,ii) - c62*s2(:,ii) + c63*s3(:,ii) + c64*s4(:,ii) - c65*s5(:,ii) ) );
        %w( :,ii+1 ) = w( :,ii ) + h*( (25/216)*s1(:,ii) + (1408/2565)*s3(:,ii) + (2197/4104)*s4(:,ii) - (1/5)*s5(:,ii) );
        z( :,ii+1 ) = z(:,ii) + h*( cz1*s1(:,ii) + cz2*s3(:,ii) + cz3*s4(:,ii) - cz4*s5(:,ii) + cz5*s6(:,ii) ) ;
        s7(:,ii) = fun( t(ii) + h , z(:,ii+1) ) ;
        e( :,ii+1 ) = h * abs( ce1*s1(:,ii) - ce2*s3(:,ii) + ce3*s4(:,ii) - ce4*s5(:,ii) + ce5*s6(:,ii) - ce6*s7(:,ii) ) ;
        rErr = e( :,ii + 1 ) / abs( z( :,ii + 1 ) ) ;
        if rErr <= rTol
            t( ii+1 ) = t ( ii ) + h ;
            ii = ii + 1 ;
            h = hIN ;
        else
            h = h * .8 * ( ( rTol * abs( z(:,ii+1) ) ) ./ e(:,ii+1) ).^(1/6) ;
        end
        
        y = z ;
        jj = jj + 1 ;
    end
end