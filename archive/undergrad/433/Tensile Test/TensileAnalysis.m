clear ; close all ; clc ;

load( 'RawData.mat' ) ;
A = 12.6*.96 ;
L = 150 ;
figure
hold on
for ii = 8:8
    SampleData = RawData{ii} ;
    time = SampleData( :,1 ) ;
    extension = SampleData( :,2 ) ;
    load = SampleData( :,3 ) ;
    strain = [ 0 ; SampleData( :,4 ) ] ;
    SampleData(:,2) = SampleData(:,2) - SampleData(1,2) ;
    indexB = find( SampleData( : , 2 ) <= .30 & SampleData( : , 2 ) >= .29 ) ;
    m = ( SampleData(indexB,3) - SampleData(1,3) )/SampleData(indexB,2) ;
    preExt = SampleData(1,3)/m ;
    correctedData = SampleData(:,2) + preExt ;
    correctedData = [ 0 , 0 ; correctedData , SampleData(:,3) ] ;
%     plot( SampleData(:,2) , SampleData(:,3) )
%     plot( correctedData(:,1) , correctedData(:,2) )
    [ MaxForce , IndexMaxForce ] = max( correctedData(:,2) ) ;
    fprintf( 'The max force for sample %i is %f N at %f mm of extension \n' , ii , MaxForce , correctedData(IndexMaxForce,1) ) 
    stress = correctedData(:,2)/A ;
    strainTheo = correctedData(:,1)/L ;
    plot( strainTheo , stress )
    [ MaxStress , IndexMaxStress ] = max( stress ) ;
    fprintf( 'The max stress for sample %i is %f MPa at %f percent strain \n' , ii , MaxForce , correctedData(IndexMaxForce,1) ) 
% plot( time , strain )
    if strain( floor( end/8 ) ) > 0 
        [ MaxStrain , IndexMaxStrain ] = max( strain ) ;
        plot( strain(1:IndexMaxStrain) , stress(1:IndexMaxStrain) )
        plot( strainTheo , stress )
    end
    if strain( floor( end/8 ) ) < 0 
        [ MinStrain , IndexMinStrain ] = min( strain ) ;
        plot( strain(1:IndexMinStrain) , stress(1:IndexMinStrain) )
    end
end
    xlabel( 'Strain' )
    ylabel( 'strain [MPa]' )
%     xlabel( 'Extension [mm]' )
%     ylabel( 'Load [N]' )
    axis( [ -.01 , .01 , 0 , 300 ] )
hold off