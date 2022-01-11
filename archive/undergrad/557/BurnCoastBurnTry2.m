clear ; close all ; clc ;
T = .1 ;
ve = .9 ;   
x = 1.05 ;
y = 0 ; 
xdot = 0 ;
ydot = sqrt( 1/x ) ;
mass = 1 ;
lambda = [ -2 ; -1 ; -1 ; -2 ; -1 ] ;
lambda = [-0.546630073119638,-0.147919753839651,-0.159151280132354,-0.771166552733118,-0.708676128788879]' ;
sto0 = [ 2 ; 0 ; 0 ; 1/sqrt(2) ; 1 ; lambda ] ;
tg = 6 ;
tol = 1e-10 ;
guess = [ lambda ; tg ] ;
opts = optimoptions( 'fsolve' , 'Display' , 'iter-detailed'  , 'FunctionTolerance' , tol ,  ...
    'OptimalityTolerance' , tol , 'MaxIterations' , 1e6 , 'MaxFunctionEvaluations' , 1e6 , 'Algorithm' , 'levenberg-marquardt' ) ; % 
optsode = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[ guess , F ] = fsolve( @BCBfpFun , guess , opts ) ;
s0fp = [ x ; y ; xdot ; ydot ; mass ; guess(1:5) ] ;
[ tfp , sfp ] = ode45( @BCBEOM , [ 0 , guess(6) ] , s0fp , optsode , T , ve ) ;
    [ tso , sso ] = ode45( @BCBEOM , [ 0 , 10 ] , s0fp , optsode , 0 , ve ) ;
    [ teo , seo ] = ode45( @BCBEOM , [ 0 , 30 ] , sfp(end,1:10) , optsode , 0 , ve ) ;
lamfp = sfp(:,6:10) ;
for ii = 1:length( tfp )
    lamv = sqrt( lamfp(ii,3)^2 + lamfp(ii,4)^2 ) ;
    switchingfp(ii) = (lamv)*ve + lamfp(ii,5)*sfp(ii,5) ;
end
figure
hold on
plot( tfp , switchingfp )
plot( [ 0 , max(tfp) ] , [ 0 , 0 ] , '-r' )
hold off
 
figure
axis equal
hold on
plot( sso(:,1) , sso(:,2) )
plot( sfp(:,1) , sfp(:,2) )
plot( seo(:,1) , seo(:,2) )

% lambda = ' ;
lambda = [ -2 ; -3 ; -2 ; -2 ; -1 ] ;
sto0 = [ 2 ; 0 ; 0 ; 1/sqrt(2) ; 1 ; lambda ] ;
t1 = 2.6 ;
t2 = 4.7 - t1 ; 
tf = 4.7*1.2 - t2 - t1 ;
guess = [ lambda ; t1 ; t2 ; tf ] ;
[ guess , F ] = fsolve( @BCBSolveFun , guess , opts ) ;
s0 = [ x ; y ; xdot ; ydot ; mass ; guess(1:5) ] ;
t1 = abs( guess(6) ) ;
t2 = abs( guess(7) ) + t1 ;
tf = abs( guess(8) ) + t2 ;
[ tb1 , sb1 ] = ode45( @BCBEOM , [ 0 , t1 ] , s0(1:10) , opts , T , ve ) ;
[ tc , sc ] = ode45( @BCBEOM , [ t1 , t2 ] , sb1(end,:) , opts , 0 , ve ) ;
[ tb2 , sb2 ] = ode45( @BCBEOM , [ t2 , tf ] , sc(end,:) , opts , T , ve ) ;
t = [ tb1 ; tc ; tb2 ] ;
s = [ sb1 ; sc ; sb2 ] ;
lam = s(:,6:10) ;
    [ tso , sso ] = ode45( @BCBEOM , [ 0 , 10 ] , s0(1:10) , optsode , 0 , ve ) ;
    [ teo , seo ] = ode45( @BCBEOM , [ 0 , 30 ] , s(end,1:10) , optsode , 0 , ve ) ;
for ii = 1:length( t )
    lamv = sqrt( lam(ii,3)^2 + lam(ii,4)^2 ) ;
    switching(ii) = lamv*ve + lam(ii,5)*s(ii,5) ;
end
figure
hold on
plot( t , switching )
plot( [ 0 , max(t) ] , [ 0 , 0 ] , '-r' )
hold off
fprintf( 'The mass percentage used for constant thrust is %f \n' , 1 - sfp(end,5) )
fprintf( 'The mass percentage used for burn-coast-burn is %f \n' , 1 - s(end,5) ) 
figure
axis equal
hold on
plot( s(:,1) , s(:,2) ) 
plot( sso(:,1) , sso(:,2) )
plot( seo(:,1) , seo(:,2) )
legend( 'Transfer' , 'Starting Orbit' , 'Final Orbit' )


function F = BCBfpFun( guess )
    T = .1 ;
    ve = .9 ;
    lam0 = guess(1:5) ;
    tf = abs( guess(end) ) ;
        x = 1.05 ;
        y = 0 ; 
        xdot = 0 ;
        ydot = sqrt( 1/x ) ;
        mass = 1 ;
    s0 = [ x ; y ; xdot ; ydot ; mass ] ;
    
    optsode = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [ t , s ] = ode45( @BCBEOM , [ 0 , tf ] , [ s0 ; lam0 ] , optsode , T , ve ) ;

    lamf = s(end,6:10) ;
    sf = s(end,1:5) ;
    
    
    % omega 
    rf = sqrt( sf(1)^2 + sf(2)^2 ) ;
%     F(1,1) = ( sf(1)^2 + sf(2)^2 - 4 )*100 ;
%     F(2,1) = ( sf(3)^2 + sf(4)^2 - 1/2 )*100 ;
%     F(3,1) = ( sf(4)*sf(1) - sf(3)*sf(2) - 2/sqrt(2) )*100 ;
    F(1,1) = ( sf(1)^2 + sf(2)^2 - 4 ) ;
    F(2,1) = ( sf(3)^2 + sf(4)^2 - 1/2 ) ;
    F(3,1) = ( sf(4)*sf(1) - sf(3)*sf(2) - 2/sqrt(2) ) ;
    
    % lambda final
    F(4,1) = lamf(5) + 1 ;
    
    % Switching
        lam0v = sqrt( lam0(3)^2 + lam0(4)^2 ) ;
    F(5,1) = ((lam0v)*ve + lam0(5)*s0(5)) ;
        lamfv = sqrt( lamf(3)^2 + lamf(4)^2 ) ;
    F(6,1) = ( (lamfv)*ve + lamf(5)*sf(5) ) ;
    
    % dependency
    F(7,1) = -lamf(1)*( sf(2)/sf(1) ) + lamf(2) - lamf(3)*( sf(4)/sf(3) )*...
        ( ( sf(3) + ( sf(2)/sf(1) )*sf(4) )/( sf(1) + ( sf(4)/sf(3) )*sf(2) ) )...
        + lamf(4)*( ( sf(3) + ( sf(2)/sf(1) )*sf(4) )/( sf(1) + ( sf(4)/sf(3) )*sf(2) ) ) ;
    
    % Hf
    lamv = sqrt( lamf(3)^2 + lamf(4)^2 ) ;
    c1 = -lamf(3)/lamv ;
    c2 = -lamf(4)/lamv ;
%     F(8,1) = lamf(1)*sf(3) + lamf(2)*sf(4) + lamf(3)*( -sf(1)/rf^3 )*c1 ...
%         + lamf(4)*( -sf(2)/rf^3 )*c2 ;
    F(8,1) = lamf(1)*sf(3) + lamf(2)*sf(4) + lamf(3)*( -sf(1)/rf^3 ) ...
        + lamf(4)*( -sf(2)/rf^3 ) ;
end

function F = BCBSolveFun( guess )
    T = .1 ;
    ve = .9 ;
    t = guess(6:8) ;
    t1 = abs( t(1) ) ;
    t2 = abs( t(2) ) + t1 ; 
    tf = abs( t(3) ) + t2 ;
    lam0 = guess( 1:5 ) ;
        x = 1.05 ;
        y = 0 ; 
        xdot = 0 ;
        ydot = sqrt( 1/x ) ;
        mass = 1 ;
    s0 = [ x ; y ; xdot ; ydot ; mass ] ;
    
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [ tb1 , sb1 ] = ode45( @BCBEOM , [ 0 , t1 ] , [ s0 ; lam0 ] , opts , T , ve ) ;
    [ tc , sc ] = ode45( @BCBEOM , [ t1 , t2 ] , sb1(end,:) , opts , 0 , ve ) ;
    [ tb2 , sb2 ] = ode45( @BCBEOM , [ t2 , tf ] , sc(end,:) , opts , T , ve ) ;
    t = [ tb1 ; tc ; tb2 ] ; 
    s = [ sb1 ; sc ; sb2 ] ;
    lamb = sb1(end,6:10) ;
    lamc = sc(end,6:10) ;
    lamf = sb2(end,6:10) ;
    sb = sb1(end,1:5) ;
    sc = sc(end,1:5) ;
    sf = sb2(end,1:5) ;
    
    % omega 
    rf = sqrt( sf(1)^2 + sf(2)^2 ) ;
    vf = sqrt( 1/rf ) ;
    F(1,1) = ( sf(1)^2 + sf(2)^2 - 4 );
    F(2,1) = ( sf(3)^2 + sf(4)^2 - 1/2  ) ;
    F(3,1) = ( sf(4)*sf(1) - sf(3)*sf(2) - 2/sqrt(2) ) ;
    
    % lambda final
    F(4,1) = lamf(5) + 1;
    
    % Switching
        lam0v = sqrt( lam0(3)^2 + lam0(4)^2 ) ;
    F(5,1) = (lam0v*ve + lam0(5)*s0(5)) ;
        lambv = sqrt( lamb(3)^2 + lamb(4)^2 ) ;
    F(6,1) = (lambv*ve + lamb(5)*sb(5)) ;
        lamcv = sqrt( lamc(3)^2 + lamc(4)^2 ) ;
    F(7,1) = (lamcv*ve + lamc(5)*sc(5)) ;
        lamfv = sqrt( lamf(3)^2 + lamf(4)^2 ) ;
    F(8,1) = (lamfv*ve + lamf(5)*sf(5)) ;
    
    % dependency
    F(9,1) = -lamf(1)*( sf(2)/sf(1) ) + lamf(2) - lamf(3)*( sf(4)/sf(3) )*...
        ( ( sf(3) + ( sf(2)/sf(1) )*sf(4) )/( sf(1) + ( sf(4)/sf(3) )*sf(2) ) )...
        + lamf(4)*( ( sf(3) + ( sf(2)/sf(1) )*sf(4) )/( sf(1) + ( sf(4)/sf(3) )*sf(2) ) ) ;
    
    % Hf
    rf = sqrt( sf(1)^2 + sf(2)^2 ) ;
    lamv = sqrt( lamf(3)^2 + lamf(4)^2 ) ;
    c1 = -lamf(3)/lamv ;
    c2 = -lamf(4)/lamv ;
    F(10,1) = lamf(1)*sf(3) + lamf(2)*sf(4) + lamf(3)*( -sf(1)/rf^3 ) ...
        + lamf(4)*( -sf(2)/rf^3 ) ;
end

function ds = BCBEOM( t , s , T , ve )
    lam = s(6:10) ;
    lamv = sqrt( lam(3)^2 + lam(4)^2 ) ;
    c1 = -lam(3)/lamv ;
    c2 = -lam(4)/lamv ;
    r = sqrt( s(1)^2 + s(2)^2 ) ;
    
    % f
    ds(1,1) = s(3) ;
    ds(2,1) = s(4) ;
    ds(3,1) = -s(1)/r^3 + (T/s(5))*c1 ;
    ds(4,1) = -s(2)/r^3 + (T/s(5))*c2 ;
    ds(5,1) = -T/ve ;
    
    % lambda dot
    ds(6,1) = lam(3)/r^3 - ( 3*s(1)/r^5 )*( lam(3)*s(1) + lam(4)*s(2) ) ;
    ds(7,1) = lam(4)/r^3 - ( 3*s(2)/r^5 )*( lam(4)*s(2) + lam(3)*s(1) ) ;
    ds(8,1) = -lam(1) ;
    ds(9,1) = -lam(2) ;
    ds(10,1) = ( T/s(5)^2 )*( lam(3)*c1 + lam(4)*c2 ) ;
end
