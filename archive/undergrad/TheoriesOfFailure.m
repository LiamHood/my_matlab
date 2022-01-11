clear ; close all ; clc ;
sigy = [ 303.4 ; 972.2 ; 1103.2 ; 302 ]*10^6 ;
signt = [ 206.46 ; 23.905 ; 1.1334 ] ;
sigmt = [ 235.45 ; 29.367 ; 8.475 ] ;
sigft = [ 254.52 ; .46 ; -.31361 ] ;
FSnt = ToF( sigy , signt ) ;
FSmt = ToF( sigy , sigmt ) ;
FSft = ToF( sigy , sigft ) ;

function FS = ToF( sigy , sigp ) 
FS = zeros( length( sigy ) , 2 ) ;
    for ii = 1:length( sigy )
        trescay = sigy(ii)/2 ;
        trescap = ( max(sigp) - min(sigp) )/2 ;
        FS(ii,1) = trescay/trescap ;
        vmy = 2*sigy(ii)^2 ;
        vmp = (sigp(1)-sigp(2))^2 + (sigp(2)-sigp(3))^2 + (sigp(3)-sigp(1))^2 ;
        FS(ii,2) = sqrt( vmy / vmp ) ;
    end
end
