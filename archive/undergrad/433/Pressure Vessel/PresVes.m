clear ; close all ; clc ;
% % Data{1} = xlsread( "Specimen_RawData_1.csv" ) ;
% % Data{2} = xlsread( "Specimen_RawData_2.csv" ) ;
% % Data{3} = xlsread( "Specimen_RawData_3.csv" ) ;
% % Data{4} = xlsread( "Specimen_RawData_4.csv" ) ;
% % Data{5} = xlsread( "Specimen_RawData_5.csv" ) ;
% % Data{6} = xlsread( "Specimen_RawData_6.csv" ) ;
% % Data{7} = xlsread( "Specimen_RawData_7.csv" ) ;
% % Data{8} = xlsread( "Specimen_RawData_8.csv" ) ;
% % Data{9} = xlsread( "Specimen_RawData_9.csv" ) ;
% % save( "Data.mat")
% 
E = 69e9 ;
t = .11e-3 ;
r = 33.605e-3 ;
v = .33 ;
Area = r^2*pi - ( r - t )^2*pi ;
% 
% load( "Data.mat" )
% for ii = 1:length( Data ) 
%     time = Data{ii}(:,1) ;
%     extension = Data{ii}(:,2) ;
%     load = Data{ii}(:,3) ;
%     strain = Data{ii}(:,4) ;
%     stress = load/Area ;
%     
%     jj = 0 ;
%     ds = 0 ;
%     while ds <= .001
%         jj = jj + 1 ;
%         ds = abs( strain(jj+1) - strain(jj) ) ;
%     end
%     istrain = average( strain( jj-50:jj-2 ) ) ;
%     fstrain = average( strain( jj+10:jj+50 ) ) ;
%     
%     figure
%     plot( time , strain )
% %     figure
% %     plot( strain , stress )

hoopstrain = [ 0.00118 , 0.00104 , 0.00147 ] ;
longstrain = [ 2.8e-04 , 3.5e-04 , 2.7e-04 , 2.7e-04 , 3.9e-04 ] ;
hstrainm = mean( hoopstrain ) ;
lstrainm = mean( longstrain ) ;
hstrains = std( hoopstrain ) ;
lstrains = std( longstrain ) ;

hstress = ( E/( 1 - v^2 ) )*( hstrainm + v*lstrainm ) ;
lstress = ( E/( 1 - v^2 ) )*( lstrainm + v*hstrainm ) ;
Ph = t*hstress/r ;
Phs = 2*t*E*hstrainm/( r*( 2 - v ) ) ;
Pl = 2*t*lstress/r ;
Pls = 2*t*E*lstrainm/( r*( 1 - 2*v ) ) ;
PercentDiffInPress = (Pl - Ph ) / Pl ;