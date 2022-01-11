function [ C ] = DavenportQ( w , sb , sg )
    
    Bt = zeros( 3,3 ) ;
    for ii = 1:length(w)
        Bt = Bt + w(ii)*sb(:,ii)*sg(:,ii)' ;
    end
    B = Bt' ;
    k22 = trace( B ) ;
    k12 = [ B(2,3) - B(3,2) ; B(3,1) - B(1,3) ; B(1,2) - B(2,1) ] ;
    k11 = B + B' - k22*eye(3) ;
    K = [ k11 , k12 ; k12' , k22 ] ;
    [ evec , lam ] = eig( K ) ;
    [ ~ , column ] = max( max( lam ) ) ;
    qbook = evec( : , column )' ;
    q = [ qbook(4) , qbook(1:3) ] ;
    C = quat2rotm( q ) ;
    disp( 'Why is DavenQ transposed? ' )
    
end