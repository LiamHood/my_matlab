function [lst] = LSidereal( time , date , location )
%Calculates Local Siderial Time 
%   Input arguments are time in [hour,minute,second] and date in
%   [day,month,year] and location is in east longitude. West longitude is
%   negative
J2000 = 2451545 ;
JC = 36525 ;
    [ ~ , Jo , UT ] = Julian( time , date ) ; %Runs Julian to find Jo
    To = ( Jo - J2000 )/JC ;
    theta_G0 = 100.4606184 + 36000.77004*To + .000387933*To^2-2.58e3*To^3 ;
    if theta_G0 > 360 % If to big to 
        while theta_G0 > 360
            theta_G0 = theta_G0 - 360 ;
        end
    elseif theta_G0 < 0
        while theta_G0 < 0
            theta_G0 = theta_G0 + 360 ;
        end
    end
    theta_G = theta_G0 + 360.98564724*(UT/24) ;
    theta = theta_G + location ;
        if theta > 360
        while theta > 360
            theta = theta - 360 ;
        end
    elseif theta < 0
        while theta < 0
            theta = theta + 360 ;
        end
    end
    lst = theta ;
    
end

