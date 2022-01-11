clear ; close all ; clc ;
load( 'RawData.mat' )

loadTheo = linspace( 0 , 5 , 1e4 ) ;
figure
hold on
for ii = 5
    data{ii} =  [ zeros( 1 , 4 ) ; data{ii} ] ;
    time = data{ii}(:,1) ;
    extension = data{ii}(:,2) ;
    load = data{ii}(:,3) ;
    sg = data{ii}(:,4) ;
    poly = polyfit( load , sg , 1 ) ;
    strainTheo = loadTheo*poly(1) ;
    sg = sg - poly(2) ;
    plot( sg , load , '.')
    plot( strainTheo , loadTheo )
end
hold off
strainAt4 = strainTheo( find( loadTheo <= 4.0005 & loadTheo > 4 ) )  
strainAt4real = sg( 254 ) 
figure
hold on
for ii = 5
    data{ii} =  [ zeros( 1 , 4 ) ; data{ii} ] ;
    time = data{ii}(:,1) ;
    extension = data{ii}(:,2) ;
    load = data{ii}(:,3) ;
    sg = data{ii}(:,4) ;
    poly = polyfit( load , extension , 1 ) ;
    extensionTheo = loadTheo*poly(1) ;
    extension = extension - poly(2) ;
    plot( extension , load , '.')
    plot( extensionTheo , loadTheo )
end
hold off
extensionAt4 = extensionTheo( find( loadTheo <= 4.0005 & loadTheo > 4 ) )  


E = 69e9 ; 
loadTheo = linspace( 0 , 5 , 1e4 ) ;
figure
hold on
for ii = 5
    data{ii} =  [ zeros( 1 , 4 ) ; data{ii} ] ;
    time = data{ii}(:,1) ;
    extension = data{ii}(:,2) ;
    load = data{ii}(:,3) ;
    sg = data{ii}(:,4) ;
    poly = polyfit( load , sg , 1 ) ;
    strainTheo = loadTheo*poly(1) ;
    sg = sg - poly(2) ;
    sigma = E*sg ;
    plot( sg , sigma )
end
hold off
stressAt4 = sigma( find( sg <= strainAt4real & sg > strainAt4real + strainAt4real*1e-3 ) )  