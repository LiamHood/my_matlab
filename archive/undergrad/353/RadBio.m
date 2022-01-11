clear ; clc ; close all ;
photons = linspace( 6 , 6 , 36 ) ; % mGy
electrons = linspace( 6 , 6 , 36 ) ; % mGy
protons = linspace( 6 , 6 , 36 ) ; % mGy
protons(8) = 20 ;
electrons(8) = 20 ;
protons(17) = 30 ;
electrons(17) = 30 ;
monthlydose = linspace( 6 , 18 , 7 ) ;
protons(22:28) = monthlydose ;
electrons(22:28) = monthlydose ;
photons(22:28) = monthlydose ;
absorbed = photons + electrons + protons ;

Hphot = photons*1 ;
Hele = electrons*1 ;
Hpro = protons*2 ;
H = Hphot+Hele+Hpro ;

CumuDoseL = zeros(1,37) ;
for ii = 1:36
   CumuDoseL(ii+1) = CumuDoseL(ii) + H(ii) ;
end
CumuDose = CumuDoseL( 2:37 ) ;

MonthDoses = [ zeros( 1 , 11 ) , H ] ;
YearlyDose = zeros( 1 , 36 ) ;
for ii = 1:36
    for jj = 1:12
        YearlyDose(ii) = YearlyDose(ii) + MonthDoses(ii+jj-1) ;
    end
end

months = 1:36 ;
figure
plot( months , photons , months , electrons , months , protons , months , absorbed , months , H , months , CumuDose , months , YearlyDose )
title( 'Radiation Dose' )
xlabel( 'Mission Month' )
ylabel( 'Dose (mGy for absorbed, mSv for effective)' )
legend( 'Photon Absorbed Dose' , 'Electrons Absorbed Dose' , 'Proton Absorbed Dose' , 'Absorbed Monthly Dose' , 'Effective Monthly Dose' , 'Total Effective Dose' , 'Yearly Effectve Dose' , 'Location' , 'northwest' ) 
disp( 'a' )
disp( 'The total effective dose only reaches 1146 mSv and the career limit is 2500 mSv.' )
disp( 'b' )
alpha = 120 ; % mGy
effectiveDose = alpha*20 ; 
disp( 'With an effective dose of 2.4 Sv there is a 35% chance of fatality in 30 days' )

