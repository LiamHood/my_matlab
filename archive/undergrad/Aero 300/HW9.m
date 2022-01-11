%% HW 9
% Liam Hood
clear
close all

%% 12.1.4
a = [ 10 -12 -6 ; 5 -5 -4 ; -1 0 3 ] ;
b = [ -14 20 10 ; -19 27 12 ; 23 -32 -13 ] ;
c = [ 8 -8 -4 ; 12 -15 -7 ; -18 26 12 ] ;
d = [ 12 -4 -2 ; 19 -19 -10 ; -35 52 27 ] ;
iteration = 10 ;
guess = [ .8 .4 .4 ]'  ;
[ veca , lama ] = rqi( a , guess , iteration ) ;
[ vecb , lamb ] = rqi( b , guess , iteration ) ;
[ vecc , lamc ] = rqi( c , guess , iteration ) ;
[ vecd , lamd ] = rqi( d , guess , iteration ) ;
guess = [ -.5 -.5 .5 ]' ;
[ veca2 , lama2 ] = rqi( a , guess , iteration ) ;
[ vecb2 , lamb2 ] = rqi( b , guess , iteration ) ;
[ vecc2 , lamc2 ] = rqi( c , guess , iteration ) ;
[ vecd2 , lamd2 ] = rqi( d , guess , iteration ) ;
guess = [ 0 .45 -.9 ]' ;
[ veca3 , lama3 ] = rqi( a , guess , iteration ) ;
[ vecb3 , lamb3 ] = rqi( b , guess , iteration ) ;
[ vecc3 , lamc3 ] = rqi( c , guess , iteration ) ;
[ vecd3 , lamd3 ] = rqi( d , guess , iteration ) ;
%% 12.2.5a
A = [ 4 3 1 ; -5 -3 0 ; 3 2 1 ] ;
[ Qbar , lam ] = qr_eigen( A , iteration ) 

%% 12.3
A2 = [ 0 1 ; 1 1.5 ] ;
[ v , lam1 ] = eig(A2) ;
s = abs(lam1) ;
for ii = 1:length(A2)
    if lam1(ii,ii) >= 1 
        u(:,ii) = v(:,ii) ;
    else
        u(:,ii) = -v(:,ii) ;
    end
end
disp( 'U' )
disp( u )
disp( 'S' )
disp( s )
disp( 'V' )
disp( v )
%% Functions
function [ vec , lam ]=rqi( A , guess , iteration )
    x = guess ;
    for j = 1:iteration
        u = x/norm(x) ;                % normalize
        lam = u'*A*u ;                 % Rayleigh quotient
        x = (A-lam*eye(size(A)))\u ;   % inverse power iteration
        x = x/norm(x) ;
    end
    vec = x/norm(x) ;
    lam = u'*A*u ;                     % Rayleigh quotient
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
