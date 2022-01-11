clear ; close all ; clc ;

syms th1 th2 th3 thh px py pz nx ny nz ox oy oz ax ay az ;
a2 = 9 ;
a3 = 9 ;
% px = 2 ;
% py = 3 ;
% pz = 5 ;
rotz = @(x) [ cos(x) , -sin(x) , 0 , 0 ;
              sin(x) , cos(x) , 0 , 0 ;
              0 , 0 , 1 , 0 ;
              0 , 0 , 0 , 1 ] ;
          
rotx = @(x) [ 1 , 0 , 0 , 0 ;
              0 , cos(x) , -sin(x) , 0 ;
              0 , sin(x) , cos(x) , 0 ;
              0 , 0 , 0 , 1 ] ;

transz = @(x) [ 1 , 0 , 0 , 0 ;
                0 , 1 , 0 , 0 ;
                0 , 0 , 1 , x ;
                0 , 0 , 0 , 1 ] ;

transx = @(x) [ 1 , 0 , 0 , x ;
                0 , 1 , 0 , 0 ;
                0 , 0 , 1 , 0 ;
                0 , 0 , 0 , 1 ] ;

% th1 = atan( -px/py ) ;
% th3 = acos(((py/cos(th1))^2+(pz)-162)/162) ;
A1 = [ -sin(th1) , 0 , -cos(th1) , 0 ; 
        cos(th1) , 0 , -sin(th1) , 0 ;
        0 , -1 , 0 , 0 ; 0 , 0 , 0 , 1 ] ;
A2 = [ sin(th2) , cos(th2) , 0 , 9*sin(th2) ;
        -cos(th2) , sin(th2) , 0 , -9*cos(th2) ;
        0 , 0 , 1 , 0 ; 0 , 0 , 0 , 1 ] ;
A3 = [ cos(th3) , -sin(th3) , 0 , 9*cos(th3) ;
        sin(th3) , cos(th3) , 0 , 9*sin(th3) ;
        0 , 0 , 1 , 0 ; 0 , 0 , 0 , 1 ] ;
Ah = [ cos(thh) , 0 , sin(thh) , 0 ;
        sin(thh) , 0 , -cos(thh) , 0 ;
        0 , 1 , 0 , 0 ; 0 , 0 , 0 , 1 ] ;

% A1 = rotz( 0 )*rotx( pi/2 ) ;
% A2 = rotz( 0 )*transx( 0 ) ;
% A3 = rotz( 0 )*transx( 0 ) ;
% Ah = rotz( 0 )*rotx( pi/2 ) ;

uTh = A1*(A2*(A3*Ah)) ;
% uTh = A3*Ah ;
matrix = [ nx , ox , ax , px ;
           ny , oy , ay , py ;
           nz , oz , az , pz ;
           0  , 0  , 0  , 1  ] ;
       
LHS = inv(A2)*inv(A1)*matrix ;
RHS = inv(A2)*inv(A1)*uTh ;

equ1 = LHS(1,4) == RHS(1,4) ;
equ2 = LHS(2,4) == RHS(2,4) ;
equ3 = LHS(3,4) == RHS(3,4) ;
disp(LHS(1,4))
disp(RHS(1,4))
disp(LHS(2,4))
disp(RHS(2,4))
disp(LHS(3,4))
disp(RHS(3,4))
% S = solve( [ equ1 , equ2 , equ3] , [ th1, th3 , th2] , 'ReturnConditions' , true ) ;

% a = axes ;
% fimplicit(equ1,[-2*pi 2*pi],'b')
% hold on
% grid on
% fimplicit(equ2,[-2*pi 2*pi],'m')
% fimplicit(equ3,[-2*pi 2*pi],'r')
% L = sym(-2*pi:pi/2:2*pi) ;
% a.XTick = double(L) ;
% a.YTick = double(L) ;
% M = arrayfun(@char, L, 'UniformOutput', false);
% a.XTickLabel = M;
% a.YTickLabel = M;
% title('Plot of System of Equations')
