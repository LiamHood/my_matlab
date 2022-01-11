function xdot = SpringDamper( t , x , nf , dr )
    
     xdot(1) = x(2) ;
     xdot(2) = -2*nf*dr*x(2)-nf^2*x(1) ;
     xdot = xdot' ;

end