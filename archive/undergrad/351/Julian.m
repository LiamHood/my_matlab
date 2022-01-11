function [ Jd , Jo , UT , J2000 ] = Julian( time , date )
%Calculates the Julian Date from a date and time
%   Uses an input of date in form [dd,mm,yyyy] and time in UT [hour,minute,second] to find Julian date. BCE years should be
%   negative

% Julian date without time
    Jo = 367*date(3) - floor(( 7*( date(3)+floor(( date(2)+9 )/12 )) )/4)+floor((275*date(2))/9) + date(1) + 1721013.5 ;
% Time 
    hour = time(1) ; %hours past noon as fraction of a day
    minute = time(2)/(60) ; %minutes as fraction of a day
    second = time(3)/(60*60) ; %seconds as fraction of a day
    UT = hour + minute + second ; % add time together 
    Jd = Jo + UT/24 ; % Full Julian date
    
% J2000 date
    J2000 = Jd - 2451545 ;
end