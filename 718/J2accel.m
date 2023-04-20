function aP = J2accel( r )
    mu = 398600 ;
    J2 = -1.08262617385222e-3 ;
    Re = 6378.1 ;
    aP = zeros(3,1) ;
    rmag = norm(r) ;
    aP(1) = -3*J2*mu*Re^2*r(1)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
    aP(2) = -3*J2*mu*Re^2*r(2)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
    aP(3) = -3*J2*mu*Re^2*r(3)/(2*rmag^5)*(3-5*r(3)^2/rmag^2) ;
end