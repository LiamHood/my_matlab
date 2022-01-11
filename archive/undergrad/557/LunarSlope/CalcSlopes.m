clear ; close all ; clc ;
latFull = -90:10:90 ;
latFull = -latFull ;
longFull = -180:10:180 ;
hsize = 2560 ;
R = 1737400 ;

m = length( latFull ) ;
n = length( longFull ) ;
TitleBeg = 'LDEM_' ;
TitleEnd = '.mat' ;
slopeFreqAll = 0 ;
for kk = 2:m 
    for LL = 2:n
        clear h ;
        strTitle = [ TitleBeg , num2str( latFull(kk) ) , '_' , num2str( longFull(LL) ) , TitleEnd ] ;
        load( strTitle ) ;
        h = double( h ) ;
        latSeg = linspace( latFull(kk-1) , latFull(kk) , hsize ) ;
        longSeg = linspace( longFull(kk-1) , longFull(kk) , hsize ) ;
        for ii = 1:hsize
            for jj = 1:hsize
                y(ii,jj) = latSeg(ii).*R ;
                x(ii,jj) = R.*longSeg(jj).*cosd( latSeg(ii) ) ;
            end
        end
        aa = 0 ;
        segSlope = 0 ;
        maxSlopeList = 0 ;
        avgSlopeList = 0 ;
        for ii = 2:(hsize-1)
            for jj = 2:(hsize-1) 
                aa = aa + 1 ;
                slope(1) = slopeFun( x(ii-1,jj-1) , x(ii,jj) , y(ii-1,jj-1) , y(ii,jj) , h(ii-1,jj-1) , h(ii,jj) ) ;
                slope(2) = slopeFun( x(ii,jj-1) , x(ii,jj) , y(ii,jj-1) , y(ii,jj) , h(ii,jj-1) , h(ii,jj) ) ;
                slope(3) = slopeFun( x(ii+1,jj-1) , x(ii,jj) , y(ii+1,jj-1) , y(ii,jj) , h(ii+1,jj-1) , h(ii,jj) ) ;
                slope(4) = slopeFun( x(ii-1,jj) , x(ii,jj) , y(ii-1,jj) , y(ii,jj) , h(ii-1,jj) , h(ii,jj) ) ;
                slope(5) = slopeFun( x(ii+1,jj) , x(ii,jj) , y(ii+1,jj) , y(ii,jj) , h(ii+1,jj) , h(ii,jj) ) ;
                slope(6) = slopeFun( x(ii-1,jj+1) , x(ii,jj) , y(ii-1,jj+1) , y(ii,jj) , h(ii-1,jj+1) , h(ii,jj) ) ;
                slope(7) = slopeFun( x(ii,jj+1) , x(ii,jj) , y(ii,jj+1) , y(ii,jj) , h(ii,jj+1) , h(ii,jj) ) ;
                slope(8) = slopeFun( x(ii+1,jj+1) , x(ii,jj) , y(ii+1,jj+1) , y(ii,jj) , h(ii+1,jj+1) , h(ii,jj) ) ;
%                 segSlope((8*(aa-1)+1):8*aa) = slope ;
%                 maxSlopeP(ii,jj) = max(slope) ;
%                 maxSlopeList(aa) = max( slope ) ;
%                 avgSlopeP(ii,jj) = mean(slope) ;
%                 avgSlopeList(aa) = avgSlopeP(ii,jj) ;
                avgSlopeList(aa) = mean(slope) ;
            end
        end
        favgSlope = floor( avgSlopeList.*(180/pi)*100 )/100 ;
        maxAvgSlope = max( favgSlope ) ;
        minAvgSlope = min( favgSlope ) ;
        posSlopes = -60:.01:60 ;
        for ii = 1:length( posSlopes )
            index = find( favgSlope == posSlopes(ii) ) ;
            slopeFreq(ii) = length( index ) ;
        end
        slopeFreqAll = slopeFreqAll + slopeFreq ;
%         rowh = linspace( 0 , size(h,1) , m ) ;
%         colh = linspace( 0 , size(h,2) , n ) ;
    end
end
save( 'AvgAbsSlopeFrequency.mat' , 'slopeFreq' )
% figure
% plot( posSlopes , slopeFreq ) 

function slope = slopeFun(xp,xo,yp,yo,hp,ho)
    slope = (hp-ho)/sqrt( (xp-xo)^2 + (yp-yo)^2 ) ;
    if slope <= 0 
        slope = -slope ;
    end
end
