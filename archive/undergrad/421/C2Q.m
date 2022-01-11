function [ eps , eta ] = C2Q( RotationMatrix )
    C = RotationMatrix ;
    eta = .5*sqrt( 1+trace(C) ) ;
    eps = zeros( 3 , 1 ) ;
    if eta ~= 0 
        eps(1) = (C(2,3)-C(3,2))/(4*eta) ;
        eps(2) = (C(3,1)-C(1,3))/(4*eta) ;
        eps(1) = (C(1,2)-C(2,1))/(4*eta) ;
    else
        eps(1) = sqrt( (C(1,1)+1)/2 ) ;
        eps(2) = sqrt( (C(2,2)+1)/2 ) ;
        eps(3) = sqrt( (C(3,3)+1)/2 ) ;
    end
end