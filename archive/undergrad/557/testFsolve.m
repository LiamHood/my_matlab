clear ; close all ; clc ;
% x = 0 ;
% xdot = 0 ; 
% lambda = [ -1 ; -1 ] ;
% s0 = [ x ; xdot ; lambda ] ;
% 
% opts = optimoptions( 'fsolve' , 'Display' , 'iter' , 'Algorithm' , 'levenberg-marquardt' ) ;
% [ s , F ] = fsolve( @ex1 , s0 , opts ) ;




clear ;
x = 0 ;
y = 0 ;
xdot = 0 ; 
ydot = 0 ;
lambda = [ -1 ; -1 ; -1 ; -1 ] ;
s0 = [ x ; y ; xdot ; ydot ; lambda ] ;
ttrans = 3 ;

opts = optimoptions( 'fsolve' , 'Display' , 'iter-detailed'  , 'FunctionTolerance' , 1e-8 , 'StepTolerance' , 1e-8 ) ; % , 'Algorithm' , 'levenberg-marquardt'
[ s , F ] = fsolve( @SteeringAngle1 , s0 , opts ) ;
optsode = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tspan = [ 0 , 3 ] ;
[ tall , sall ] = ode45( @SteeringAngle1EOM , tspan , s , optsode ) ;
for ii = 1:length( tall )
    bot = sqrt( s(7)^2 + s(8)^2 ) ;
    cosc(ii) = -sall(ii,7)/bot ;
%     if cosc > 0
%         cosc = -cosc ;
%     end
    sinc(ii) = -sall(ii,8)/bot ;
%     if sinc > 0
%         sinc = -sinc ;
%     end
    c(ii) = atan2d( sinc(ii) , cosc(ii) ) ;
end 
figure
hold on
plot( sall(:,1) , sall(:,2) )
quiver( sall(:,1) , sall(:,2) , cosc' , sinc' , 'AutoScaleFactor' , .25 )
hold off
figure
plot( sall(:,3) , sall(:,4) )
figure
plot( tall , c )

function F = SteeringAngle1( s0 )
    tspan = [ 0 , 3 ] ;
    
    
    opts = odeset( 'RelTol' , 1e-10 , 'AbsTol' , 1e-10 ) ;
    [ t , s ] = ode45( @SteeringAngle1EOM , tspan , s0 , opts ) ;
    sf = s(end,:) ;
    
    % force intitial conditions
    F(1,1) = s0(1) ;
    F(2,1) = s0(2) ;
    F(3,1) = s0(3) ;
    F(4,1) = s0(4) ;
    
    % omega 
    F(5,1) = sf(2) - 1 ;
    F(6,1) = sf(4) ;
    
    % lambda final
    F(7,1) = sf(5) ;
    F(8,1) = sf(7) + 1 ;
    % DONT NEED LAMBDAS EQUAL TO LITTLE OMEGAS
%     w1 = s0(6) ;
%     w2 = s0(6)*tspan(2) + s0(8) ;
    
end

function ds = SteeringAngle1EOM( t , s )
    bot = sqrt( s(7)^2 + s(8)^2 ) ;
    cosc = -s(7)/bot ;
    sinc = -s(8)/bot ;

    ds = zeros( 8 , 1 ) ;
    ds(1) = s(3) ;
    ds(2) = s(4) ;
    ds(3) = cosc ;
    ds(4) = sinc ;
    ds(5) = 0 ;
    ds(6) = 0 ;
    ds(7) = -s(5) ;
    ds(8) = -s(6) ;
end

function F = ex1( s0 )    
    w1 = s0(3) ;
    w2 = -.5*s0(3) ;
    tspan = [ 0 , 1 ] ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [ t , s ] = ode45( @ex1EOM , tspan , s0 , opts ) ;
    sf = s(end,:) ;
    F = zeros( length(s0)*1.5 , 1 ) ;
    % omega 
    F(1) = sf(1) - 1 ;
    F(2) = sf(2) ;
    
    % lambda final
    F(3) = sf(3) - w1 ;
    F(4) = sf(4) - w2 ;
    
    % force intitial conditions
%     F(5) = s0(1) ;
%     F(6) = s0(2) ;
%     F(7) = s0(3) ;
%     F(8) = s0(4) ;
    
end

function ds = ex1EOM( t , s )
    
    ds = zeros( 4 , 1 ) ;
    ds(1) = s(2) ;
    ds(2) = -s(4) ;
    ds(3) = 0 ;
    ds(4) = -s(3) ;
   
end