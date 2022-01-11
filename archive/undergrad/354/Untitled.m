opts = odeset( 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 );
s = 9.74 ;
m = 10 ;
A = (2*2.54*1e-2)^2 ;
tspan = [ 0 , 10 ] ;
T0 = 293*ones(3,1) ;
Tsource = 500 ;
[ t , Tnew ] = ode45( @MLI , tspan , T0 , opts , eps , A , Tsource , s , m , abs ) ;
figure
plot( t , Tnew(:,1) , t , Tnew(:,2) , t , Tnew(:,3)  )
function [dT,t] = MLI( t , T , eps , A , Tsource , s , m , abs )
sig = 5.6704e-8 ;
qdotto1 = abs*Tsource^4*A ;
dT(1) = qdotto1/(s*m) ;
    for ii = 1:( length( T ) - 1 )
        qdot(ii) = (eps*abs*sig*(T(ii)-T(ii+1))^4)/(2*((1-eps)/(A*eps))+(1/A)) ;
    end
    dT(2:length(T)) = (qdot(1:(length(T)-1))./(s*m)) ;
    dT = dT' ;
end