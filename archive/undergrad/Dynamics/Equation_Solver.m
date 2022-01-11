%example using the Matlab symbolic solver to solve a system of 4
%simultaneous equations
%Birdsong 2/13/14
syms abx aby aab abc %define the variable names for symbolic solution
%% define the system of equations and solve using the 'solve' command
% note you don't need to arrange them in any particular order unlike 
% using linear algebra and inverting a matrix

S = solve( abx == -3 * ( -0.8787 )^2 - 3.9 * aab , ...
    aby == -3.9 * ( -0.8787 )^2 + 3 * aab , ...
    0 == abx - 2.3 * 1.146^2 + 0.5 * abc , ...
    0 == aby + 0.5 * 1.146^2 + 2.3 * abc ) ;

S=[S.abx, S.aby, S.aab, S.abc]; %extract the answers from the S object
%vpa(S,4) %convert answer to a 4 digit number

%% Change this section
syms w_ab w_bc %define the variable names for symbolic solution

A = solve( 0 == w_ab^2 * 3 + w_bc^2 * 2.3 , ...
    4 == w_ab *  - w_bc * 2.3 );

A=[ A.w_ab , A.w_bc ] %extract the answers from the S object
vpa(A,2) %convert answer to a 4 digit number
