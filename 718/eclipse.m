function [T, ecl_time, ecl_proportion] = eclipse(a, mu, Rb)
% Finds orbital period, eclipse time in seconds, and eclipse proportion of
% orbit
% Needs inputs of semi-major axis, mu of central body, radius of central
% body

    beta = asin(Rb/a);
    T = 2*pi*sqrt(a^3/mu);
    ecl_proportion = 2*beta/(2*pi);
    ecl_time = ecl_proportion*T;
    ecl_min = ecl_time/60;
    
    fprintf("The period is %f minutes\n", T/60)
    fprintf("The proportion of eclipse is %f \n", ecl_proportion)
    fprintf("The eclipse is %f minutes\n", ecl_min)
end