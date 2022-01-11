clear
close all

%fun = @(t,y) ( y - t - 1 )^2 + 2 ;
fun = @(t,y) y^2 + ( -2 - 2*t ) * y + 3 + 2*t + t^2 ;
tspan = [ 0 pi/3 ] ;
y0 = 0 ;
h = .01 ;
rTol = 10^-6 ;

[t,w,z,e,eother,rErr] = rkf45(fun, tspan, y0, h, rTol) ;

[ tODE , yODE ] = ode45( fun , tspan , y0 );
figure

yAct = @(t) tan(t) + t + 1 ;
gErrODE = abs( yODE - yAct( tODE ) ) ;
plot( t , w , '-b' , t , z , '-r' )
figure
plot( t , eother , '-r' )
figure
plot( t , e , '-b' )

function [t,y] = rkf45(fun, tspan, y0, h, rTol)
% Runge-Katta 4 with 5th order check
w(:,1) = y0 ;
z(:,1) = y0 ;
t(1) = tspan(1) ;
ii = 1 ;
jj = 1 ;
hIN = h ;
    while t(ii) < tspan(2)
        s1(:,ii) = fun( t(ii) , w(:,ii) ) ;
        s2(:,ii) = fun( t(ii) + (1/4)*h , w(:,ii) + (1/4)*h*s1(:,ii) ) ;
        s3(:,ii) = fun( t(ii) + (3/8)*h , w(:,ii) + (3/32)*h*s1(:,ii) + (9/32)*s2(:,ii) ) ;
        s4(:,ii) = fun( t(ii) + (12/13)*h , w(:,ii) + (1932/2197)*h*s1(:,ii) - (7200/2197)*h*s2(:,ii) + (7296/2197)*h*s3(:,ii) ) ;
        s5(:,ii) = fun( t(ii) + h , w(:,ii) + (439/216)*h*s1(:,ii) - (8)*h*s2(:,ii) + (3680/513)*h*s3(:,ii) - (845/4104)*h*s4(:,ii)) ;
        s6(:,ii) = fun( t(ii) + .5*h , w(:,ii) - (8/27)*h*s1(:,ii) + (2)*h*s2(:,ii) - (3544/2565)*h*s3(:,ii) + (1859/4104)*h*s4(:,ii) - (11/40)*h*s5(:,ii) ) ;
        w( :,ii+1 ) = w( :,ii ) + h*( (25/216)*s1(:,ii) + (1408/2565)*s3(:,ii) + (2197/4104)*s4(:,ii) - (1/5)*s5(:,ii) );
        z( :,ii+1 ) = z(:,ii) + h*( (16/135)*s1(:,ii) + (6656/12825)*s3(:,ii) + (28561/56430)*s4(:,ii) - (9/50)*s5(:,ii) + (2/55)*s6(:,ii) ) ;
        e( :,ii+1 ) = h * abs( (1/360)*s1(:,ii) - (128/4275)*s2(:,ii) - (2197/75240)*s4(:,ii) + (1/50)*s5(:,ii) + (2/55)*s6(:,ii) ) ;
        eother( :,ii+1 ) = abs( w( :,ii+1 ) - z( :,ii+1 ) );
        rErr = e( :,ii + 1 ) / abs( w( :,ii + 1 ) ) ;
        if rErr <= rTol
            t( ii+1 ) = t ( ii ) + h ;
            w( :,ii+1 ) = z( :,ii+1 ) ;
            ii = ii + 1 ;
            h = hIN ;
        else
            h = h * .5 * ( ( rTol*abs(w(ii+1)) ) / e(ii+1) )^(1/5) ;
        end
        
        y = w ;
        jj = jj + 1 ;
    end
end

% function [t,w,z,e,eother,rErr] = rkf45(fun, tspan, y0, h, rTol)
% % Runge-Katta 4 with 5th order check
% w(1) = y0 ;
% z(1) = y0 ;
% t(1) = tspan(1) ;
% ii = 1 ;
% jj = 1 ;
% hIN = h ;
%     while t(ii) < tspan(2)
%         s1 = fun( t(ii) , w(ii) ) ;
%         s2 = fun( t(ii) + (1/4)*h , w(ii) + (1/4)*h*s1 ) ;
%         s3 = fun( t(ii) + (3/8)*h , w(ii) + (3/32)*h*s1 + (9/32)*s2 ) ;
%         s4 = fun( t(ii) + (12/13)*h , w(ii) + (1932/2197)*h*s1 - (7200/2197)*h*s2 + (7296/2197)*h*s3 ) ;
%         s5 = fun( t(ii) + h , w(ii) + (439/216)*h*s1 - (8)*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4) ;
%         s6 = fun( t(ii) + .5*h , w(ii) - (8/27)*h*s1 + (2)*h*s2 - (3544/2565)*h*s3 + (1859/4104)*h*s4 - (11/40)*h*s5 ) ;
%         w( ii+1 ) = w( ii ) + h*( (25/216)*s1 + (1408/2565)*s3 + (2197/4104)*s4 - (1/5)*s5 );
%         z( ii+1 ) = z(ii) + h*( (16/135)*s1 + (6656/12825)*s3 + (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6 ) ;
%         e( ii+1 ) = h * abs( (1/360)*s1 - (128/4275)*s2 - (2197/75240)*s4 + (1/50)*s5 + (2/55)*s6 ) ;
%         eother( ii+1 ) = abs( w( ii+1 ) - z( ii+1 ) );
%         rErr = e( ii + 1 ) / abs( w( ii + 1 ) ) ;
%     
%             t( ii+1 ) = t ( ii ) + h ;
%             ii = ii + 1 ;
%             h = hIN ;
%       
%         
%         jj = jj + 1 ;
%     end
% end
