function i = sso(a, e)
    mu = 398600;
    J2 = -1.08262617385222e-3 ;
    Re = 6378.1;
    
    T = sqrt(a^3/mu);
    T_ES = 365.249*24*3600;
    cos_i = 2*pi*(T/T_ES)*a^2*(1-e^2)^2/(-3*pi*J2*Re^2);
    i = acos(cos_i);
end