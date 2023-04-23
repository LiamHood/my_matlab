function i = sso(a, e)
    mu = 398600;
    J2 = -1.08262617385222e-3 ;
    Re = 6378.1;
    dRAAN_sunsyn = 1.991063853e-7;
    
    cos_i = (-2*a^(7/2)*dRAAN_sunsyn*(1-e^2)^2)/(3*Re^2*J2*sqrt(mu));
    i = acos(cos_i);
end