function [ pf ] = f2d( p )
%Uses a given vector to create a polynomial of any size then finds all
%non-zero vectors
    lp = length(p); %saves the length of the input vector
    pf = [ p ; zeros( lp - 1 , lp )  ] ; %allocates a space for the derivatives below the original vector
    fm = zeros( lp ) ; %allocates space for the function matrix

    %Finds the coeffecient of each derivative. Each row down is another
    %derivative
    for ii = 2:lp %only adding information to the second row and below
        p = polyder( p ) ; %Takes derivative of last row
        for jj = 0:( length(p) - 1 ) %fills the values of p into the main matrix from right to left 
        pf( ii , lp - jj ) = p( length(p) - jj ) ;
        end
    end

%     %Uses the coeffecients from the matrix this evaluates each row for a value
%     %of x
%     for kk = 1:lp %
%         for mm = 1:lp
%         fm(kk,mm) = ( x^( lp - mm ) * pf( kk , mm ) ); %Evaluates each cell by multiplying the value by the appropriate value of x
%         end
%     end
%     f = sum(fm') ; %adds all the terms of the function matrix
end

