%% Aero 215 HW2
% Liam Hood
% 10/13/17
% Question 4


%% Calculations
mach_h = 20.05 * (288.16) ^ .5 * 3.28084 ; %speed of sound at 10,000 ft in ft/s
mach_l = 20.05 * (268.348) ^ .5 * 3.28084; %speed of sound sea level in ft/s
z = .5 * ( mach_h * .5 )^2 / 32.174 + 10000 ; %finding energy height
disp( 'Energy height of mach 0.5 at 10,000 ft' )
disp( z )
v = ( ( z - 10000 ) * 2 * 32.174 ) ^ .5 ; %Finding the speed at sea level for the above enrgy height
v_m = v / mach_l ; %Turns speed into a mach number
disp( 'Mach' )
disp( v_m )
disp( 'Speed of sound at 10,000 ft in ft/s' )
disp( mach_h )
disp( 'Speed of sound at sea level in ft/s' )
disp( mach_l )
