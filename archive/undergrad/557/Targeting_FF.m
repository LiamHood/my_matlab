function [ va , vb , deltav , collision ] = Targeting_FF( rint , rtgt , vint , vtgt , dt , mu , Rcb )
% finds low deltav lamberts solution from a fixed point state to a fixed
% state. No propagation of target through transit time. Checks for
% collision with central body
    Tran_n = cross( rint , rtgt ) ;
    h_n = cross( rint , vint ) ;
    if dot( Tran_n , h_n ) > 0
        tm = 1 ;
    else 
        tm = -1 ;
    end
    [ va , vb ] = Lambert_Izzo( rint , rtgt , dt , tm , mu ) ;
    collision = HitCB( rint , rtgt , va , vb , mu , Rcb) ;
    deltav = abs( norm( va - vint ) + norm( vb - vtgt ) ) ;
    
    function [ collision ] = HitCB( rint , rtgt , va , vb , mu , Rcb )
        if dot( rint , va ) > 0 && dot( rtgt , vb )
            coes = state2coes( [ rint ; va ] , mu ) ;
            if coes(8) < Rcb
                collision = 1 ;
            else
                collision = 0 ;
            end
        else
            collision = 0 ;
        end
    end
end