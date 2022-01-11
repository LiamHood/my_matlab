function ap = SRPacceleration( r , v , A , m , JDo , t , solarsail) 
% if solar sail is one then it increases the area by 100 when sun and
% velocity are aligned
    ii = 0 ;
    p = 4.57e-6 ;
    [ light , s ] = ValladoShadow( r , JDo , t ) ;
    if solarsail == 1 
        alignment = ( dot( v , s )/( norm(s)*norm(v) ) ) ;
        if alignment < -.7
            A = A*100 ;
        end
    end 
    cr = 1.2 ;
    ap = (- ( p * cr * A * s * light )/m)*1e-3 ;
    
end