clear ; close all ; clc ;
% bottles1 = xlsread( 'Lab6g1' ) ;
% bottles2 = xlsread( 'Lab6g2' ) ;
% bottles3 = xlsread( 'Lab6g3' ) ;
% bottles4 = xlsread( 'Lab6g4' ) ;
% save( 'Bottles.mat' )
load( 'Bottles.mat' )

for ii = 1:floor( length(bottles3)/60 ) 
    bottles3m( ii , : ) = bottles3( 60*(ii-1)+1 , : ) ;
end
for ii = 1:floor( .6*length(bottles4)/6 ) 
    bottles4m( ii , : ) = bottles4( 6*(ii-1)+2 , : ) ;
end

time{1} = bottles1( : , 1 ) ;
time{2} = bottles2( : , 1 ) ;
time{3} = bottles3m( : , 1 ) ;
time{4} = bottles4m( : , 1 ) ;

silver{1} = bottles1( : , 3 ) ;
silver{2} = bottles2( : , 3 ) ;
silver{3} = bottles3m( : , 2 ) ;
silver{4} = bottles4m( : , 3 ) ;


white{1} = bottles1( : , 2 ) ;
white{2} = bottles2( : , 4 ) ;
white{3} = bottles3m( : , 3 ) ;
white{4} = bottles4m( : , 2 ) ;

black{1} = bottles1( : , 4 ) ;
black{2} = bottles2( : , 2 ) ;
black{3} = bottles3m( : , 4 ) ;
black{4} = bottles4m( : , 4 ) ;

ambient{1} = bottles1( : , 5 ) ;
ambient{2} = bottles2( : , 5 ) ;
ambient{3} = bottles3m( : , 5 ) ;
ambient{4} = bottles4m( : , 5 ) ;

sig = 5.67037441e-8 ;
r = 2*2.54 ;
h = 8*2.54 ;
area = pi*r*h ; % surface area assuming half lit at all times
temperr = .5 ; % degrees C of Error
volerr = 10 ;
s = 4.186 ; % J/gram*K
m = 300 ; % grams

for ii = 1:4 
    err = ones( length(time{ii}) , 1 )*temperr ;
    errqs = zeros( length(time{ii}) , 1 ) ;
    errqw = zeros( length(time{ii}) , 1 ) ;
    errqb = zeros( length(time{ii}) , 1 ) ;
    qSilver = zeros( length(time{ii}) , 1 ) ;
    qWhite = zeros( length(time{ii}) , 1 ) ;
    qBlack = zeros( length(time{ii}) , 1 ) ;
    for jj = 1:length( time{ii} )
        qSilver(jj) = s*m*( silver{ii}(jj) - min( silver{ii} ) ) ;
            errqs(jj,1) = ( 2*temperr./silver{ii}(jj) + (volerr/m) ).*qSilver(jj) ;
        qWhite(jj) = s*m*( white{ii}(jj) - min( white{ii} ) ) ;
            errqw(jj,1) = ( 2*temperr./white{ii}(jj) + (volerr/m) ).*qWhite(jj) ;
        qBlack(jj) = s*m*( black{ii}(jj) - min( black{ii} ) ) ;
            errqb(jj,1) = ( 2*temperr./black{ii}(jj) + (volerr/m) ).*qBlack(jj)   ;
    end
    figure
    hold on
%     errorbar( time{ii} , silver{ii} , err , 'b' )
%     errorbar( time{ii} , white{ii} , err , 'r' )
%     errorbar( time{ii} , black{ii} , err , 'k' )
%     errorbar( time{ii} , ambient{ii} , err )
%     legend( 'Bare Aluminum' , 'White' , 'Black' ,'Ambient' , 'Location' , 'southeast' )
%     title( 'Surface Treatment Effect on Water Heat' )
%     hold off
%     
    errorbar( time{ii} , qSilver , errqs , 'b' )
    errorbar( time{ii} , qWhite , errqw , 'r' )
    errorbar( time{ii} , qBlack , errqb , 'k' )
    legend( 'Bare Aluminum' , 'White' , 'Black' , 'Location' , 'southeast' )
    title( 'Surface Treatment Effect on Energy In' )
    xlabel( 'Time (s)' )
    ylabel( 'Thermal Energy In (J)' )
    hold off
end