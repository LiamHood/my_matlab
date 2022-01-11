clear ; close all ; clc ;
% LDEM = Tiff('Lunar_LDEM_118m.tif','r') ;
% imshow(imageData) ;

% imageData = read(LDEM) ;
% load(  'SampleDEM.mat' )
load(  'LDEM_118m_NW.mat' )
% imageDataO = imageData ;
% imageData = imageData( 1:100 , 1:200 ) ;
R = 1737400 ;
mpy = 118.4505876 ;
[ m , n ] = size( imageData ) ;
latall = linspace( -pi/2 , 0 , m ) ;
longall = linspace( 0 , pi , n ) ;
for ii = 1:n
    for jj = 1:m
        y(jj,ii) = latall(jj).*R ;
        x(jj,ii) = R.*longall(ii).*cos( latall(jj) ) ;
    end
end
h = imageData.*0.5 ;
figure
surf( x , y , h )

