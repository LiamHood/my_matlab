function [ S , C ] = Stumpf( z ) 
% Stumpff Functions
sl = 11 ;
sc = zeros(1,sl) ;
selement = zeros(1,sl) ;
cc = zeros(1,sl) ;
celement = zeros(1,sl) ;
S = 0 ;
C = 0 ;

        for kk = 1:sl
            sc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+3) ;
        end
        for kk = 1:sl
            cc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+2) ;
        end
        
        for jj = 1:sl
            selement(jj) = sc(jj)*z^(jj-1) ;
        end
        S = sum(selement) ;
        
        for jj = 1:sl
            celement(jj) = cc(jj)*z^(jj-1) ;
        end
        C = sum(celement) ;
end