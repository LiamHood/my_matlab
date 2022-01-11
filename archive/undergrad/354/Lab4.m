clear ; close all ; clc ;
% data = xlsread( 'Lab4_354.xlsx' ) ;
load( 'Lab4_Group1' ) ;
load( 'Lab4_Group2' ) ;
load( 'Lab4_Group3' ) ;
load( 'Lab4_Group4' ) ;
Tsun = data( 3,1 ) ;
Tcloud = data( 3,7 ) ;
psun = 1366*(.06*.0375) ;

    
    isno{1} = data( 3:48 , 2 ) ;
    vsno{1} = data( 3:48 , 3 )*10^(-3) ;
    isclean{1} = data( 3:41 , 5 ) ;
    vsclean{1} = data( 3:41 , 6 )*10^(-3) ;
    icno{1} = data( 3:24 , 8 ) ;
    vcno{1} = data( 3:24 , 9 )*10^(-3) ;
    vcno{1}(6,1) = .4830 ;
    icdirty{1} = data( 3:19 , 10 ) ;
    vcdirty{1} = data( 3:19 , 11 )*10^(-3) ;
    icpoly{1} = data( 3:21 , 12 ) ;
    vcpoly{1} = data( 3:21 , 13 )*10^(-3) ;
    
    isno{2} = [ data2( 1:16 , 1 ) ; data2( 1:14 , 7 ) ; data2( 1:14 , 13 ) ]  ;
    vsno{2} = [ data2( 1:16 , 2 ) ; data2( 1:14 , 8 ) ; data2( 1:14 , 14 ) ]  ;
    isclean{2} = [ data2( 25:38 , 1 ) ; data2( 25:37 , 7 ) ; data2( 25:36 , 13 ) ]  ;
    vsclean{2} = [ data2( 25:38 , 2 ) ; data2( 25:37 , 8 ) ; data( 25:36 , 14 ) ]  ;
    icno{2} = [ data2( 125:140 , 1 ) ; data2( 125:139 , 7 ) ; data2( 125:140 , 13 ) ]  ;
    vcno{2} = [ data2( 125:140 , 2 ) ; data2( 125:139 , 8 ) ; data2( 125:140 , 14 ) ]  ;
    icdirty{2} = [ data2( 102:114 , 1 ) ; data2( 102:113 , 7 ) ; data2( 102:87 , 13 ) ]  ;
    vcdirty{2} = [ data2( 102:114 , 2 ) ; data2( 102:113 , 8 ) ; data2( 102:87 , 14 ) ]  ;
    icpoly{2} = [ data2( 77:93 , 1 ) ; data2( 77:95 , 7 ) ; data2( 77:96 , 13 ) ]  ;
    vcpoly{2} = [ data2( 77:93 , 2 ) ; data2( 77:95 , 8 ) ; data2( 77:96 , 14 ) ]  ;
    
    isno{3} = data3( 1:34 , 2 ) ;
    vsno{3} = data3( 1:34 , 1 )*10^(-3) ;
    isclean{3} = data3( 1:34 , 6 ) ;
    vsclean{3} = data3( 1:34 , 5 )*10^(-3) ;
    icno{3} = [ data2( 125:140 , 1 ) ; data2( 125:139 , 7 ) ; data2( 125:140 , 13 ) ]  ;
    vcno{3} = [ data2( 125:140 , 2 ) ; data2( 125:139 , 8 ) ; data2( 125:140 , 14 ) ]  ;
    icdirty{3} = data3( 1:34 , 10 ) ;
    vcdirty{3} = data3( 1:34 , 9 )*10^(-3) ;
    icpoly{3} = data3( 1:30 , 14 ) ;
    vcpoly{3} = data3( 1:30 , 13 )*10^(-3) ;
    
    isno{4} = data( 3:48 , 2 ) ;
    vsno{4} = data( 3:48 , 3 )*10^(-3) ;
    isclean{4} = data4( 1:36 , 1 ) ;
    vsclean{4} = data4( 1:36 , 2 )*10^(-3) ;
    icno{4} = data( 3:24 , 8 ) ;
    vcno{4} = data( 3:24 , 9 )*10^(-3) ;
    icdirty{4} = data4( 1:35 , 7 ) ;
    vcdirty{4} = data4( 1:35 , 8 )*10^(-3) ;
    icpoly{4} = data4( 1:21 , 13 ) ;
    vcpoly{4} = data4( 1:21 , 14 )*10^(-3) ;
    
for ii = 1:4 
    isnosc{ii} = max(isno{ii}) ;
    vsnooc{ii} = max(vsno{ii}) ;
    iscleansc{ii} = max( isclean{ii} ) ;
    vscleanoc{ii} = max( vsclean{ii} ) ;
    icnosc{ii} = max( icno{ii} ) ;
    vcnooc{ii} = max( vcno{ii} ) ;
    icdirtysc{ii} = max( icdirty{ii} ) ;
    vcdirtyoc{ii} = max( vcdirty{ii} ) ;
    icpolysc{ii} = max( icpoly{ii} ) ;
    vcpolyoc{ii} = max( vcpoly{ii} ) ;
    ptsno{ii} = isnosc{ii} * vsnooc{ii} ;
    ptsclean{ii} = iscleansc{ii} * vscleanoc{ii} ;
    ptcno{ii} = icnosc{ii} * vcnooc{ii} ;
    ptcdirty{ii} = icdirtysc{ii} * vcdirtyoc{ii} ;
    ptcpoly{ii} = icpolysc{ii} * vcpolyoc{ii} ;
    errorI = .01 ; % smallest measurement was .01 Amps
    errorV = .005 ; % voltages varied that changes of up to 5 volts could go unnoticed    

    [ prsno{ii} , ind ] = max( isno{ii}.*vsno{ii} ) ;
    eprsno{ii} = (errorI/isno{ii}(ind) + errorV/vsno{ii}(ind))*prsno{ii} ;
    [ prsclean{ii} , ind ] = max( isclean{ii}.*vsclean{ii} ) ;
    eprsclean{ii} = (errorI/isclean{ii}(ind) + errorV/vsclean{ii}(ind))*prsclean{ii} ;
    [ prcno{ii} , ind ] = max( icno{ii}.*vcno{ii} ) ;
    eprcno{ii} = (errorI/icno{ii}(ind) + errorV/vcno{ii}(ind))*prcno{ii} ;
    [ prcdirty{ii} , ind ] = max( icdirty{ii}.*vcdirty{ii} ) ;
    eprcdirty{ii} = (errorI/icdirty{ii}(ind) + errorV/vcdirty{ii}(ind))*prcdirty{ii} ;
    [ prcpoly{ii} , ind ] = max( icpoly{ii}.*vcpoly{ii} ) ;
    eprcpoly{ii} = (errorI/icpoly{ii}(ind) + errorV/vcpoly{ii}(ind))*prcpoly{ii} ;

    
    ffsno{ii} = prsno{ii}./ptsno{ii} ;
        effsno{ii} = ( eprsno{ii}/prsno{ii} + errorI/isnosc{ii} + errorV/vsnooc{ii} )*ffsno{ii} ;
    ffsclean{ii} = prsclean{ii}./ptsclean{ii} ;
        effsclean{ii} = ( eprsclean{ii}/prsclean{ii} + errorI/iscleansc{ii} + errorV/vscleanoc{ii} )*ffsclean{ii} ;
    ffcno{ii} = prcno{ii}./ptcno{ii} ;
        effcno{ii} = ( eprcno{ii}/prcno{ii} + errorI/icnosc{ii} + errorV/vcnooc{ii} )*ffcno{ii} ;
    ffcdirty{ii} = prcdirty{ii}./ptcdirty{ii} ;
        effcdirty{ii} = ( eprcdirty{ii}/prcdirty{ii} + errorI/icdirtysc{ii} + errorV/vcdirtyoc{ii} )*ffcdirty{ii} ;
    ffcpoly{ii} = prcpoly{ii}./ptcpoly{ii} ;
        effcpoly{ii} = ( eprcpoly{ii}/prcpoly{ii} + errorI/icpolysc{ii} + errorV/vcpolyoc{ii} )*ffcpoly{ii} ;
    
    eArea = ( .001/.06 + .001/.0375 ) ;
    etasno{ii} = prsno{ii}./psun ;
        eetasno{ii} = ( eprsno{ii}/prsno{ii} + eArea )*etasno{ii} ;
    etasclean{ii} = prsclean{ii}./psun ;
        eetasclean{ii} = ( eprsclean{ii}/prsclean{ii} + eArea )*etasclean{ii} ;
    etacno{ii} = prcno{ii}./psun ;
        eetacno{ii} = ( eprcno{ii}/prcno{ii} + eArea )*etacno{ii} ;
    etacdirty{ii} = prcdirty{ii}./psun ; 
        eetacdirty{ii} = ( eprcdirty{ii}/prcdirty{ii} + eArea )*etacdirty{ii} ;
    etacpoly{ii} = prcpoly{ii}./psun ;
        eetacpoly{ii} = ( eprcpoly{ii}/prcpoly{ii} + eArea )*etacpoly{ii} ;
end

% figure 
% hold on 
% for ii = 1:4
% 
% plot( ii , etasno{ii} , '*b' ) 
% plot( ii , etasclean{ii} , '*r' ) 
% plot( ii , etacno{ii} , '*y') 
% plot( ii , etacdirty{ii} , '*g' ) 
% plot( ii , etacpoly{ii} , '*k' ) 
% 
% end
% hold off

% figure 
% hold on 
% for ii = 1:4
% 
% plot( ii , ffsno{ii} , '*b' ) 
% plot( ii , ffsclean{ii} , '*r' ) 
% plot( ii , ffcno{ii} , '*y') 
% plot( ii , ffcdirty{ii} , '*g' ) 
% plot( ii , ffcpoly{ii} , '*k' ) 
% 
% end
% hold off
ffsnosum = 0 ;
for ii = 1:4
    ffsnosum = ffsno{ii} + ffsnosum ;
end
ffsnoavg = ffsnosum/4 ;
etasnosum = 0 ;
for ii = 1:4
    etasnosum = etasno{ii} + etasnosum ;
end
etasnoavg = etasnosum/4 ;
ffscleansum = 0 ;
for ii = 1:4
    ffscleansum = ffsclean{ii} + ffscleansum ;
end
ffscleanavg = ffscleansum/4 ;
etascleansum = 0 ;
for ii = 1:4
    etascleansum = etasclean{ii} + etascleansum ;
end
etascleanavg = etascleansum/4 ;
ffcnosum = 0 ;
for ii = 1:4
    ffcnosum = ffcno{ii} + ffcnosum ;
end
ffcnoavg = ffcnosum/4 ;
etacnosum = 0 ;
for ii = 1:4
    etacnosum = etacno{ii} + etacnosum ;
end
etacnoavg = etacnosum/4 ;
% ffcdirtysum = 0 ;
% for ii = 1:4
%     ffcdirtysum = ffcdirty{ii} + ffcdirtysum ;
% end
% ffcdirtyavg = ffcdirtysum/4 ;
% etacdirtysum = 0 ;
% for ii = 1:4
%     etacdirtysum = etacdirty{ii} + etacdirtysum ;
% end
% etacdirtyavg = etacdirtysum/4 ;
% ffcpolysum = 0 ;
% for ii = 1:4
%     ffcpolysum = ffcpoly{ii} + ffcpolysum ;
% end
% ffcpolyavg = ffcpolysum/4 ;
% etacpolysum = 0 ;
% for ii = 1:4
%     etacpolysum = etacpoly{ii} + etacpolysum ;
% end
% etacpolyavg = etacpolysum/4 ;
figure

hold on
    errorIlong = errorI*ones(1,length(icno{1})) ; 
    errorVlong = errorV*ones(1,length(vcno{1})) ; 
errorbar( vcno{1} , icno{1} , errorIlong , errorIlong , errorVlong , errorVlong , 'b.' )
    errorIlong = errorI*ones(1,length(icdirty{1})) ; 
    errorVlong = errorV*ones(1,length(vcdirty{1})) ;
errorbar( vcdirty{1} , icdirty{1} , errorIlong , errorIlong , errorVlong , errorVlong , 'r.' )
    errorIlong = errorI*ones(1,length(icpoly{1})) ; 
    errorVlong = errorV*ones(1,length(vcpoly{1})) ;
errorbar( vcpoly{1} , icpoly{1} , errorIlong , errorIlong , errorVlong , errorVlong , 'g.' )
    errorIlong = errorI*ones(1,length(icpoly{2})) ; 
    errorVlong = errorV*ones(1,length(vcpoly{2})) ;
errorbar( vcpoly{2} , icpoly{2} , errorIlong , errorIlong , errorVlong , errorVlong , 'g.' )
    errorIlong = errorI*ones(1,length(icpoly{4})) ; 
    errorVlong = errorV*ones(1,length(vcpoly{4})) ;
errorbar( vcpoly{4} , icpoly{4} , errorIlong , errorIlong , errorVlong , errorVlong , 'g.' )
    errorIlong = errorI*ones(1,length(icno{2})) ; 
    errorVlong = errorV*ones(1,length(vcno{2})) ; 
errorbar( vcno{2} , icno{2} , errorIlong , errorIlong , errorVlong , errorVlong , 'b.' )
legend( 'No glass' , 'Outgassed' , 'Polymerized' , 'Location' , 'southwest' )
title( 'Cloudy Solar Cells' )
xlabel( 'Potential (V)' )
ylabel( 'Current (A)' )
hold off

figure
hold on
    errorIlong = errorI*ones(1,length(isno{1})) ; 
    errorVlong = errorV*ones(1,length(vsno{1})) ; 
errorbar( vsno{1} , isno{1} , errorIlong , errorIlong , errorVlong , errorVlong , 'b.' )
    errorIlong = errorI*ones(1,length(isclean{1})) ; 
    errorVlong = errorV*ones(1,length(vsclean{1})) ; 
errorbar( vsclean{1} , isclean{1} , errorIlong , errorIlong , errorVlong , errorVlong , 'r.' )
    errorIlong = errorI*ones(1,length(icdirty{2})) ; 
    errorVlong = errorV*ones(1,length(vcdirty{2})) ; 
errorbar( vcdirty{2} , icdirty{2} , errorIlong , errorIlong , errorVlong , errorVlong , 'g.' )
    errorIlong = errorI*ones(1,length(icpoly{3})) ; 
    errorVlong = errorV*ones(1,length(vcpoly{3})) ; 
errorbar( vcpoly{3} , icpoly{3} , errorIlong , errorIlong , errorVlong , errorVlong , 'k.' )
for ii = 3
    errorIlong = errorI*ones(1,length(icdirty{ii})) ; 
    errorVlong = errorV*ones(1,length(vcdirty{ii})) ; 
errorbar( vcdirty{ii} , icdirty{ii} , errorIlong , errorIlong , errorVlong , errorVlong , 'g.' )
end
for ii = 1:3
    errorIlong = errorI*ones(1,length(isclean{ii})) ; 
    errorVlong = errorV*ones(1,length(vsclean{ii})) ; 
errorbar( vsclean{ii} , isclean{ii} , errorIlong , errorIlong , errorVlong , errorVlong , 'r.' )    
end
for ii = 1:4
    errorIlong = errorI*ones(1,length(isno{ii})) ; 
    errorVlong = errorV*ones(1,length(vsno{ii})) ; 
errorbar( vsno{ii} , isno{ii} , errorIlong , errorIlong , errorVlong , errorVlong , 'b.' )    
end  
    errorIlong = errorI*ones(1,length(icpoly{3})) ; 
    errorVlong = errorV*ones(1,length(vcpoly{3})) ; 
errorbar( vcpoly{3} , icpoly{3} , errorIlong , errorIlong , errorVlong , errorVlong , 'k.' )
legend( 'No glass' , 'Clean Glass' ,'Outgassed' , 'Polymerized' , 'Location' , 'southwest' )
title( 'Sunny Solar Cells' )
xlabel( 'Potential (V)' )
ylabel( 'Current (A)' )
hold off
deta = zeros(1,4) ;
deta(1) = etacno{1} - etacdirty{1} ;
for ii = 2:4
    deta(ii) = etasno{ii} - etacdirty{ii} ;
end
deta(4) = etacno{1} - etacdirty{4} ;
mean( deta )
