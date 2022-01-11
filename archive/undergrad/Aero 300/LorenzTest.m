clear

param = [ 10 , 28 , 8/3 ] ;
fun = @lorenz ;
[ t , pos ] = dp45( fun , [ 0 50 ] , [ 1,1,1 ] , .01 , 10^-5 ) ;
[ tM , posM ] = ode45(fun,[0 50 ],[1,0,10]);
%plot3( pos(1,:) , pos(2,:) , pos(3,:) )
figure
plot3( posM(:,1) , posM(:,2) , posM(:,3) )


function [tt,y] = rkf45(fun, tspan, y0, h, rTol)
%Runge-Katta_Fehlberg order 4/5
    t = tspan(1) ;
    tt(:,1) = t ;
    w(:,1) = y0 ;
    ii = 1 ;
    hIn = h ;
    while t+h < tspan(2)
        s1 = fun( t , w(:,ii) ) ;
        s2 = fun( t + (1/4)*h , w(:,ii) + (1/4)*h.*s1 ) ;
        s3 = fun( t + (3/8)*h , w(:,ii) + (3/32)*h.*s1 + (9/32)*h.*s2 ) ;
        s4 = fun( t + (12/13)*h , w(:,ii) + (1932/2197)*h.*s1 - (7200/2197)*h.*s2 + (7296/2197)*h.*s3 ) ;
        s5 = fun( t + h , w(:,ii) + (439/216)*h.*s1 - 8*h.*s2 + (3680/513)*h.*s3 - (845/4104)*h.*s4 ) ;
        s6 = fun( t + .5*h , w(:,ii) - (8/27)*h.*s1 + 2*s2 - (3544/2565)*h.*s3 + (1859/4104)*h.*s4 - (11/40)*h.*s5 ) ;
        w(:,ii+1) = w(:,ii) + h.*( (25/216)*s1 + (1408/2565)*s3 + (2197/4104)*s4 - .2*s5 ) ;
        z = w(:,ii) + h.*( (16/135)*s1 + (6656/12825)*s3 + (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6 ) ;
        e = abs( z - w(:,ii+1) ) ;
        re(ii) = e/abs(w(ii+1)) ;
        if re <= rTol
            w(:,ii+1) = z ;
            t = t + h ;
            tt(:,ii+1) = t ;
            h = hIn ;
            ii = ii+1 ;
            
        else
            h = h .* .5 * ((rTol*norm(w(:,ii+1)))/e).^.2 ;
        end
    end
    y = w ;
end

function [t,y] = dp45(fun, tspan, y0, h, rTol)
%DP 4 with 5th order check
w(:,1) = y0 ;
z(:,1) = y0 ;
for ii = 1:length(y0)
    t(ii,1) = tspan(1) ;
end
ii = 1 ;
jj = 1 ;
hIN = h ;
    while t(ii) < tspan(2)
        s1(:,ii) = fun( t(:,ii) , z(:,ii) ) ;
        s2(:,ii) = fun( t(:,ii) + (1/5)*h , z(:,ii) + (1/5)*h.*s1(:,ii) ) ;
        s3(:,ii) = fun( t(:,ii) + (3/10)*h , z(:,ii) + (3/40)*h.*s1(:,ii) + (9/40).*s2(:,ii) ) ;
        s4(:,ii) = fun( t(:,ii) + (4/5)*h , z(:,ii) + (44/45)*h.*s1(:,ii) - (56/15)*h.*s2(:,ii) + (32/9)*h.*s3(:,ii) ) ;
        s5(:,ii) = fun( t(:,ii) + (8/9)*h , z(:,ii) + h.*( (19372/6561).*s1(:,ii) - (25360/2187).*s2(:,ii) + (64448/6561).*s3(:,ii) - (212/729).*s4(:,ii) ) ) ;
        s6(:,ii) = fun( t(:,ii) + h , z(:,ii) + h.*( (9017/3168).*s1(:,ii) - (355/33).*s2(:,ii) + (46732/5247).*s3(:,ii) + (49/176).*s4(:,ii) - (5103/18656).*s5(:,ii) ) );
        %w( :,ii+1 ) = w( :,ii ) + h*( (25/216)*s1(:,ii) + (1408/2565)*s3(:,ii) + (2197/4104)*s4(:,ii) - (1/5)*s5(:,ii) );
        z( :,ii+1 ) = z(:,ii) + h.*( (35/384).*s1(:,ii) + (500/1113).*s3(:,ii) + (125/192).*s4(:,ii) - (2187/6784).*s5(:,ii) + (11/84).*s6(:,ii) ) ;
        s7(:,ii) = fun( t(:,ii) + h , z(:,ii+1) ) ;
        e( :,ii+1 ) = h .* abs( (71/57600)*s1(:,ii) - (71/16695)*s3(:,ii) + (71/1920)*s4(:,ii) - (17253/339200)*s5(:,ii) + (22/525)*s6(:,ii) - (1/40)*s7(:,ii) ) ;
        rErr = norm( e( :,ii + 1 )) ./ norm( z( :,ii + 1 ) ) ;
        if rErr <= rTol
            t( :,ii+1 ) = t ( :,ii ) + h ;
            ii = ii + 1 ;
            h = hIN ;
        else
            h = h .* .8 .* ( ( rTol .* abs( z(:,ii+1)) ) ./ e(:,ii+1) ).^(1/5) ;
        end
        
        y = z ;
        jj = jj + 1 ;
    end
end

function [ pos ] = lorenz( t,pos )
%Creates lorenz equations. t is a number, pos is a row vector of [ x,y,z]
%values, param is [ sigma , rho , beta ] constants for lorenz
param = [ 10 , 28 , 8/3 ] ;
pos = [ param(1) * ( pos(2) - pos(1) ) ; pos(1) * ( param(2) - pos(3) ) - pos(2) ; pos(1) * pos(2) - param(3) * pos(3) ] ;

end