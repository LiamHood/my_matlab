function [ dstate ] = NBody( t , state , m ) 
G = 6.67259e-20 ;
L = length( state ) ;
dstate = zeros( L , 1 ) ;
t
% body state
for ii = 1:(L/6)
    r{ii} = state( (6*ii-5):(6*ii-3) ) ;
    v{ii} = state( (6*ii-2):6*ii ) ;
end

% derivative of positions is velocity
for ii = 1:(L/6)
    dstate((6*ii-5):(6*ii-3)) = v{ii};
end

% accelerations
for ii = 1:(L/6)
    for jj = 1:(L/6)
        if ii ~= jj
            dist = (norm( r{jj}-r{ii} ))^3 ;
            dstate(6*ii-2) = dstate(6*ii-2) + G*m(jj)*( r{jj}(1)-r{ii}(1) )/dist ;
            dstate(6*ii-1) = dstate(6*ii-1) + G*m(jj)*( r{jj}(2)-r{ii}(2) )/dist ;
            dstate(6*ii) = dstate(6*ii) + G*m(jj)*( r{jj}(3)-r{ii}(3) )/dist ;
        else
            dstate(6*ii-2) = dstate(6*ii-2) ;
            dstate(6*ii-1) = dstate(6*ii-1) ;
            dstate(6*ii) = dstate(6*ii) ;
        end
    end
end


end