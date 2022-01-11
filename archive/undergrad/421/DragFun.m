function Tatmo = DragFun( Cbg , vG , rG , Areas , Positions , Norms )
if norm( vG ) >= 100 
    rG = rG/1e3 ;
    vG = vG/1e3 ;
end
% Velocity and Position are pulled in as an ECI value

vnorm = norm(vG)*1000 ; % m/s Magnitude of velocity

alt = (norm(rG) - 6378) ; %km - altitude of spacecraft

% The normal vectors on the body frame are pulled in

% Use rotation matrix to rotate from ECI to Bodyframe

vb = Cbg*vG ; %rotates the velocity vector into the body frame

% Model of Atmosphere

rho = EarthStdAtm(alt) ; %inputs the altitude of the sc to get density of air

% Find the amount of area on each side, exposed to the sun, use dot products

Drag = 0 ; %preallocates force summation

Top = [0;0;0] ; %Preallocates Cps Parts 

Bot = 0 ; 

for pp = 1:length(Areas)

    dotted = (Norms(:,pp)'*(vb/norm(vb))) ; %Dot product of norm and velocity

    if dotted > 0 

        Drag = Drag + -rho*vnorm^2*dotted*Areas(pp)*(vb/norm(vb)) ;

    end 

    % Need Center of Pressure Now

    if dotted > 0 

        Top = Top + Positions(:,pp) * dotted * Areas(pp) ; % Sums the position vector times the dotted value and area

        Bot = Bot + dotted * Areas(pp) ; %Sums the n dot s * area

    end

end

Cps = Top/Bot ; %Defines the center of pressure

Tatmo = cross(Cps,Drag) ; %Cross product of force and distance

end