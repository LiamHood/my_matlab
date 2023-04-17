function aP = J3accel(r)
    mu = 398600 ;
    J3 = 2.53241051856772e-6 ;
    Re = 6378 ;
    aP = zeros(3,1) ;
    rmag = norm(r) ;
    aP(1) = -5*J3*mu*Re^3*r(1)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
    aP(2) = -5*J3*mu*Re^3*r(2)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
    aP(3) = -5*J3*mu*Re^3/(2*rmag^7)*(6*r(3)^2-7*r(3)^4/rmag^2 - (3/5)*rmag^2) ;
end