function [ Tgg ] = GravGradFun( Cbg , Ix , Iy , Iz , r )

% Include all rotations to make torque into body using an inputted State
% Finds gravity gradient torque in earth orbit

    mu = 398600*(1e3)^3 ;

    rb = (Cbg*r) ; % position vector in body frame

    rbcross = [ 0 -rb(3) rb(2) ; rb(3) 0 -rb(1) ; -rb(2) rb(1) 0 ] ;

    Ip = [ Ix , 0 , 0 ; 0 , Iy , 0 ; 0 , 0 , Iz ] ;

    Tgg = (( 3*mu )/( norm(rb)^5 ))*rbcross*Ip*rb ; % torque in body frame

end