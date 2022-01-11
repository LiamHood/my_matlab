%% Lab 10
% Liam Hood
% Aero 300

clear
close all

%% 1
%Define givens
m = 1 ;
k = .5 ;
c = .3 ;
y0 = 2 ;
v0 = 0 ;

ydot = [ 0 1 ; -k/m -c/m ] ; %Spring description
[ x , lam ] = eig(ydot) ; %Find eigen
%c1 = (v0-((x(2,2)*y0)/x(1,2)))/(x(2,1)-((x(2,2)*x(1,1))/x(1,2))) ; %Find constants using initial values and fact that when t=0 exp will always be one
%c2 = (y0-x(1,1)*c1)/x(1,2) ; 
    cons1 = inv(x)*[y0;v0] ;
    c1 = cons1(1) ;
    c2 = cons1(2) ;
%y = @(t) x(1,1)*c1*exp( real(lam(1,1))*t ).*exp( imag(lam(1,1))*t ) + x(2,1)*c2*exp( real(lam(2,2))*t ).*exp( imag(lam(2,2))*t ) ; %Displacement function
y = @(t) x(1,1)*c1*exp( ( (lam(1,1))*t ) ) + x(2,1)*c2*exp( (lam(2,2))*t ) ;
t = linspace( 0 , 50 , 1000 ) ; 
displacement = y(t) ;

%Plot of displacement
plot( t , y(t) )
title( 'Spring Displacement' )
xlabel( 'Time' )
ylabel( 'Displacement' )

%Compare to no dampening
c = 0 ;
ydot = [ 0 1 ; -k/m -c/m ] ; %Spring description
[ x2 , lam2 ] = eig(ydot) ; %Find eigen
cons = inv( x2 )*[ y0 ; v0 ] ; %Find constants using initial values and fact that when t=0 exp will always be one
y2 = @(t) x2(1,1)*cons(1)*exp( lam2(1,1)*t ) + x2(2,1)*cons(2)*exp( lam2(2,2)*t ) ; %Displacement function
t = linspace( 0 , 50 , 1000 ) ; 

%Plot of displacement with no dampening
figure
plot( t , y2(t) )
title( 'Spring Displacement with no Dampening' )
xlabel( 'Time' )
ylabel( 'Displacement' )

disp( 'If c=0 then the sin wave continues indefinitely. The eigenvalues when ' )
disp( ' c = .3 have real components and when c=0 they do not.' )
disp( 'c=.3 eigenvalues' )
disp( lam )
disp( 'c=0 eigenvalues' )
disp( lam2 )

%% 2 

% x(1) = intitial guess
% while error < tolerance
    % x(ii+1) = A*x(ii)
    % x(ii+1) = x(ii+1)/max(x(ii+1))
    % m0 = Transpose( x(ii) )*x(ii)
    % m1 = Transpose( x(ii) )*x(ii+1)
    % m2 = Transpose( x(ii+1) )*x(ii+1)
    % error = square root( (m2/m1)-(m1/m0)^2 )
% end

%% 3 and 4
%Given
A = [ 1 1 4 -4 ; -1 4 0 4 ; 4 -1 3 4 ; 1 0 2 1 ] ;

Tol = 10^-10 ; % Tolerance
[ vec3 , lam3 ] = eig(A) ; %Find eigenvalues and eigenvectors using built in function
[ vecP3 , iteration ] = eigenpower( A , Tol ) ; %Using my function
%Answer questions
disp( 'My power method function gives the largest eigenvector but with the ' )
disp( 'opposite sign as eig function.' )
disp( 'Eigenvectors followed by eigenvalues' )
disp( vec3 )
disp( lam3 )
disp( 'Eigenvector from my power method function' )
disp( vecP3 )

%% 5
%Given
A = [ 1 1 4 -4 ; -1 4 0 4 ; 4 -1 3 4 ; 1 0 2 1 ] ;
iterations = 20;
[ vecQR , lamQR ] = qr_eigen( A , iterations ) ; %Use my QR iteration function

%% 6
S = svd( A ) ; %Find singular values
scale_condition = max(S)/min(S) ; %Find scaling condition
disp( 'The scaling is relatively low so it is well conditioned' )

%% Functions
function [ vec , iteration ] = eigenpower( A , Tol )
    x = ones( length(A) , 1 ) ; %Set initial guess at all ones
    error = Tol+1 ; %Set initial error above tolerance
    iteration = 0 ; %Set iteration count at 0
    while error > Tol %Continue as long as error is above tolerance
        y = A*x ; %Calculate next guess
        y = y/norm(y) ; %Reduce size norm is one
        %Calculate values for finding error
        m0 = x'*x ; 
        m1 = x'*y ;
        m2 = y'*y ;
        error = sqrt( (m2/m1) - (m1/m0)^2 ) ; %Find error
        x = y ; %reassign newly established guess as old guess
        iteration = iteration+1 ;
    end
    vec = x ;
end

function [ Qbar , lam ] = qr_eigen( A , iterations )
%Input A of any size and iteration count. Outputs eigenvectors in Qbar and
%corresponding eigenvalues in lam
    [ height , width ] = size(A) ;
    Q = eye( height , height ) ; %Create initial guess of Q
    Qbar = Q ; %Begin Qbar
    R = A ; %Starting R value
    for jj = 1:iterations %Continues for a given number of iterations
        [ Q , R ] = qr( R*Q ) ; %Find new Q and R
        Qbar = Qbar*Q ; %Find new Qbar
    end
    lam = diag( R*Q ) ; %find eigenvalues

end

