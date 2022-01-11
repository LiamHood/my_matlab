function [ L , D ] = MarsGliderFunction( CL , CD , S , v , h )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [D0] = MarsStdAtm(h);
    q = .5 * D0 * v^2 ;
    L = CL * q * S ;
    D = CD * q * S ;
    disp( 'Lift' )
    disp( L )
    disp( 'Drag' )
    disp( D )
    disp( 'L/D' )
    disp( L/D )
   


end

