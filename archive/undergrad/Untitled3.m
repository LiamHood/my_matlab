clear ; close all ; clc ;

syms c1 c2 c3 s1 s2 s3 thh px py pz nx ny nz ox oy oz ax ay az ;
a2 = 9 ;
a3 = 9 ;
th1 = atan( -px/py ) ;
A1 = [ -s1 , 0 , -c1 , 0 ; 
        c1 , 0 , -s1 , 0 ;
        0 , -1 , 0 , 0 ; 0 , 0 , 0 , 1 ] ;
A2 = [ s2 , c2 , 0 , 9*s2 ;
        -c2 , s2 , 0 , -9*c2 ;
        0 , 0 , 1 , 0 ; 0 , 0 , 0 , 1 ] ;
A3 = [ c3 , -s3 , 0 , 9*c3 ;
        s3 , c3 , 0 , 9*s3 ;
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
       
LHS = matrix ;
RHS = uTh ;

% equx = uTh(1,4) == px ;
% equy = uTh(2,4) == py ;
% equz = uTh(3,4) == pz ;
% % a = axes ;
% % fimplicit(equx,[-2*pi 2*pi],'b')
% % hold on
% % grid on
% % fimplicit(equy,[-2*pi 2*pi],'m')
% % fimplicit(equz,[-2*pi 2*pi],'r')
% % L = sym(-2*pi:pi/2:2*pi) ;
% % a.XTick = double(L) ;
% % a.YTick = double(L) ;
% % M = arrayfun(@char, L, 'UniformOutput', false);
% % a.XTickLabel = M;
% % a.YTickLabel = M;
% % title('Plot of System of Equations')
% % 
% % % syms x y
% % % eqn1 = sin(x)+cos(y) == 4/5;
% % % eqn2 = sin(x)*cos(y) == 1/10;
% % % a = axes;
% % % fimplicit(eqn1,[-2*pi 2*pi],'b');
% % hold on
% % grid on
% % fimplicit(eqn2,[-2*pi 2*pi],'m');
% % L = sym(-2*pi:pi/2:2*pi);
% % a.XTick = double(L);
% % a.YTick = double(L);
% % M = arrayfun(@char, L, 'UniformOutput', false);
% % a.XTickLabel = M;
% % a.YTickLabel = M;
% % title('Plot of System of Equations')
% % legend('sin(x)+cos(y) == 4/5','sin(x)*cos(y) == 1/10',...
% %     'Location','best','AutoUpdate','off')
% 
% 
% S = solve( [ equx , equy , equz ] , [ th1 , th2 , th3 ] ) ;

% clear ; close all ; clc ;
% 
% syms b c x y z ;
equ1 = LHS(1,4) == RHS(1,4) ;
equ2 = LHS(2,4) == RHS(2,4) ;
equ3 = LHS(3,4) == RHS(3,4) ;
equ4 = cos(asin(s1))==c1 ;
equ5 = cos(asin(s2))==c2 ;
equ6 = cos(asin(s3))==c3 ;

S = solve( [ equ1, equ2, equ3, equ4, equ5, equ6] , [ c1, c2, c3, s1, s2, s3] ) ;
disp(S.c2(1,1))
disp(S.s2(1,1))
disp(S.s2(1,1)/S.c2(1,1))
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
