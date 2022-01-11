clear ; close all ; clc ;

data = xlsread( 'Specimen_RawData_5.csv' ) ;
time = data( :,1 ) ;
extension = data( :,2 ) ;
load = data( :,3 ) ;
sg = data( :,4 ) ;

eoff = -extension( 100 )*load(1)/(load(1)+load(100)) ;
extc = extension - eoff ;
findex = find( max( load ) == load ) ;

time = time( 1:findex ) ;
extc = extc( 1:findex ) ;
load = load( 1:findex ) ;
sg = sg( 1:findex ) ;

figure
plot( extc , load )


sfindex = find( max( sg ) == sg ) ;
sgoff = -sg( 20 )*load(1)/(load(1)+load(20)) ;
sgbf = sg( 1:sfindex ) - sgoff ;
loadsg = load( 1:sfindex ) ;
figure
plot( sgbf , loadsg )

t = 1.5 ;
w = 24 ;
D = 8.8 ; 
wd = w - D ;
A = t*w ;
Ad = t*wd ;
L0 = 128.8 ;
L0 = 140 ;
stressTheo = load/A ;
stressD = load/Ad ;
stress_sg = loadsg/A ;
stress_sgD = loadsg/Ad ;
straintheo = extc / L0 ;
figure
plot( straintheo , stressTheo , straintheo , stressD , sgbf , stress_sg , sgbf , stress_sgD ) 
legend( 'Theoretical Undamaged' , 'Theoretical Damaged' , 'Measured w/undamaged Area' , 'Measured w/damage' )