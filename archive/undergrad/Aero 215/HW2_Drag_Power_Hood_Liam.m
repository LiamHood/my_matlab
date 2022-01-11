function [ v , v_sl, D_sl, P_sl ] = HW2_Drag_Power_Hood_Liam( v_max, inputs )
%Calculates drag force, power required, and stall velocity of an aircraft
%in SLF

    v_sl = ((2 .* inputs.WS .* inputs.S ) .* ( inputs.cl_max .* inputs.rho .* inputs.S )).^.5;
    
    v = linspace( v_sl , v_max , 10 );

    D_sl = ( .5 .* inputs.Cdp .* inputs.rho .* (v.^2) .* inputs.S ) + ( inputs.W ) ./ ( .5 .* inputs.rho .* (v.^2) .* pi .* inputs.e .* inputs.AR );
    
    P_sl = ( D_sl .* v ) * .001356 ;
    
    disp( 'Velocity in ft/s' )
    disp( v )
    disp( 'Drag force in lbf' )
    disp( D_sl )
    disp( 'Power required in kW' )
    disp( P_sl )
    

end

