function C = Triad( s1a , s1b , s2a , s2b )
% s1b , s2b are vectors in body frame; s1a and s2a are vectors in inertial
% space
    xta = s1a ;
    yta = cross( s1a , s2a )/( norm( cross( s1a,s2a ))) ;
    zta = cross( xta , yta ) ;
    xtb = s1b ;
    ytb = cross( s1b , s2b )/( norm( cross( s1b,s2b ))) ;
    ztb = cross( xtb , ytb ) ;
    Cat = [ xta , yta , zta ] ;
    Cbt = [ xtb , ytb , ztb ] ;
    C = Cbt*Cat' ;
end