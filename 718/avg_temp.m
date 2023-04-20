function [T_no_ecl, T_ecl] = avg_temp(Re, Je, Js, a_orbit, a, F, F_ecl, ecl_prop, A_alpha, A_eps, alpha, eps, Q)
    sigma = 5.67e-8;
    Jp = Je*(Re/a_orbit)^2;
    Ja = Js*a*F;
    Jincident = Js + Jp + Ja;
    T_4 = A_alpha/A_eps * alpha/eps * Jincident/sigma;
    T_no_ecl = T_4^(1/4) - 272.15;
    Ja = Js*a*F_ecl;
    T_4_ecl = alpha/eps * A_alpha/A_eps * 1/sigma * [(1-ecl_prop)*(Js + Ja) + Jp + Q];
    T_ecl = T_4_ecl^(1/4) - 272.15;
    
    fprintf("No eclipse average temperature is %f degrees Celsius\n", T_no_ecl)
    fprintf("Eclipse average temperature is %f degrees Celsius\n", T_ecl)

end