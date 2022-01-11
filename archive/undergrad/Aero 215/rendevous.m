function [ delta_v ] = rendevous( r_circ , r_e , mu , ta_lead , n )
%Calculates the delta-v needed for dragon and ISS to rendevous
    T_ISS = 2 * pi * ( r_circ ^ 3 / mu ) ^ .5 ; %Period of the orbit of ISS
    T_ab = T_ISS * ( ta_lead / ( 2 * pi ) ) ; %time ISS leads dragon by
    T_phasing = T_ISS - T_ab / n ; %Period of phasing orbit
    a_phasing = ( mu * ( T_phasing / ( 2 * pi ) ) ^ 2 ) ^ ( 1 / 3 );  %semi-major axis of phasing orbit
        %Check for hitting earth
        danger = .5 * r_circ + r_e ;
        if a_phasing <= danger
            disp( 'You messed up and hit earth' )
        end
    se_phasing = - mu / ( 2 * a_phasing );  %specific energy of phasing orbit
    v_phasing = ( 2 * ( ( mu / r_circ ) + se_phasing ) ) ^ .5 ; %velocity of phasing orbit at apogee
    v_circ = ( 2 * ( mu / (2 * r_circ ) ) ) ^ .5 ; %velocity of circular ISS orbit
    delta_v = 2 * ( v_circ - v_phasing ) ; %total delta-v
end

