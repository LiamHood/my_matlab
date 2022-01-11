
clear ; close all ; clc ;

theta = 0 ;
ecc = 0 ;
a = 1.05 ;

r = a*( 1 - ecc^2 )/( 1 + ecc*cos( theta ) ) ;
rdot = ( ecc*sin( theta ) )/( sqrt( a )*sqrt( 1 - ecc^2 ) ) ;
thetadot = ( 1 + ecc*cos( theta ) )/( a*( 1 - ecc^2 ) )^1.5 ;
lambda = [ -4 ; 0 ; -1 ; -2.5 ] ;
s0 = [ r ; theta ; rdot ; thetadot ; lambda ] ;
opts = optimoptions( 'fsolve' , 'Display' , 'iter-detailed'  , 'FunctionTolerance' , 1e-8 , 'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , 'Algorithm' , 'levenberg-marquardt' ) ; % 
[ s0s , F ] = fsolve( @PolarSolveFun , s0 , opts ) ;
optsode = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[ t , s ] = ode45( @PolarEOM , [ 0 , 4 ] , s0s , optsode , .1 ) ;
[ ts , ss ] = ode45( @PolarEOM , [ 0 , 10 ] , s0s , optsode , 0 ) ;
[ te , se ] = ode45( @PolarEOM , [ 0 , 20 ] , s(end,:) , optsode , 0 ) ;
for ii = 1:length( t )
    cosc(ii) = ( -s(ii,8)/( s(ii,1)^2*s(ii,7)^2 + s(ii,8)^2 )^0.5 ) ;
    sinc(ii) = ( -( s(ii,7)*s(ii,1) )/( s(ii,1)^2*s(ii,7)^2 + s(ii,8)^2 )^0.5 ) ;
    x(ii) = s(ii,1)*cos(s(ii,2)) ;
    y(ii) = s(ii,1)*sin(s(ii,2)) ;
    c(ii) = atan2d( sinc(ii) , cosc(ii) ) ;
    if c(ii) < 0 
        c(ii) = 360 + c(ii) ;
    end
    vxmag = s(ii,3)*cos(s(ii,2)) - s(ii,1)*s(ii,4)*sin(s(ii,2)) ;
    vymag = s(ii,3)*sin(s(ii,2)) + s(ii,1)*s(ii,4)*cos(s(ii,2)) ;
    vmag = sqrt( vxmag^2 + vymag^2 ) ;
    vx(ii) = vxmag/vmag ;
    vy(ii) = vymag/vmag ;
    steer1(ii) = cosc(ii) + vx(ii) ;
    steer2(ii) = sinc(ii) + vy(ii) ;
    
end
figure
plot( t , c )
figure
polarplot( s(:,2) , s(:,1) ) 
hold on
polarplot( ss(:,2) , ss(:,1) ) 
polarplot( se(:,2) , se(:,1) ) 
hold off
figure
axis equal
hold on
plot( x , y ) 
quiver( x , y , steer1 , steer2 , 'AutoScaleFactor' , .25 )

stateLabels = [ "r" ; "theta" ; "rdot" ; "thetadot" ] ;
    stateUnits = [ "DU" ; "radians" ; "DU" ; "radians" ] ;
    for ii = 1:4
        fprintf( 'The starting %s is %f %s \n' , stateLabels(ii) , s0(ii) , stateUnits ) 
    end
    for ii = 1:4
        fprintf( 'The starting lambda %i is %f \n' , ii , s(1,ii+4) ) 
    end
    fprintf( '\n' )
    for ii = 1:4
        fprintf( 'The final %s is %f \n' , stateLabels(ii) , s(end,ii) ) 
    end
    for ii = 1:4
        fprintf( 'The final lambda %i is %f \n' , ii , s(end,ii+4) ) 
    end
    fprintf( 'The cost function is -r to maximize altitude and is %f DU \n' , -s(end,1) )
function F = PolarSolveFun( s0 )
    tspan = [ 0 , 4 ] ;
    accel = .1 ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [ t , s ] = ode45( @PolarEOM , tspan , s0 , opts , accel ) ;
    sf = s(end,:) ;
    
    % force intitial conditions
    F(1,1) = s0(1) - 1.05 ;
    F(2,1) = s0(2) ;
    F(3,1) = s0(3) ;
    F(4,1) = s0(4) - ( 1/( 1.05^1.5 ) ) ;
    
    % omega 
    F(5,1) = sf(3) ;
    F(6,1) = sf(4)^2 * sf(1)^3 - 1 ;
    
    % lambda final
    F(7,1) = sf(5) - ( -1 + 1.5*sf(8)*sf(4)/sf(1) ) ;
    F(8,1) = sf(6) ;
    F(9,1) = s0(6) ;


end
function ds = PolarEOM( t , s , accel )
    accosc = accel*( -s(8)/( s(1)^2*s(7)^2 + s(8)^2 )^0.5 ) ;
    acsinc = accel*( -( s(7)*s(1) )/( s(1)^2*s(7)^2 + s(8)^2 )^0.5 ) ;
    % f
    ds(1,1) = s(3) ;
    ds(2,1) = s(4) ;
    ds(3,1) = s(1)*s(4)^2 - ( 1/s(1)^2 ) + acsinc ;
    ds(4,1) = ( accosc - 2*s(3)*s(4) )/s(1) ;
    % lambda dot
    ds(5,1) = -s(7)*s(4)^2 - 2*s(7)/s(1)^3 + ( s(8)/s(1)^2 )*( accosc - 2*s(3)*s(4) ) ;
    ds(6,1) = 0 ;
    ds(7,1) = -s(5) + 2*s(8)*s(4)/s(1) ;
    ds(8,1) = -s(6) - 2*s(7)*s(1)*s(4) + 2*s(8)*s(3)/s(1) ;
end