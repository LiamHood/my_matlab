function [ I ] = PrincipalMoI( Ix , Iy , Iz )
I = zeros( 3,3 ) ;
I(1,1) = Ix ; 
I(2,2) = Iy ;
I(3,3) = Iz ;
end