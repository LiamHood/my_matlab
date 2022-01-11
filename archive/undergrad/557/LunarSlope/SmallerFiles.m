% clear ; close all ; clc ;
% load( 'LDEM_118m.mat' ) ;
latFull = -90:10:90 ;
latFull = -latFull ;
longFull = -180:10:180 ;
m = length( latFull ) ;
n = length( longFull ) ;
size(imageData) ;
row = linspace( 0 , size(imageData,1) , m ) ;
col = linspace( 0 , size(imageData,2) , n ) ;
TitleBeg = 'LDEM_' ;
TitleEnd = '.mat' ;
for ii = 2:m 
    for jj = 2:n 
        h = imageData( (row(ii-1)+1):row(ii) , (col(ii-1)+1):col(ii) )*.5 ;
        strTitle = [ TitleBeg , num2str( latFull(ii) ) , '_' , num2str( longFull(jj) ) , TitleEnd ] ;
        save( strTitle , 'h' )
    end
end
