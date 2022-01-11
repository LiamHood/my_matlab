function [ c ] = loc( a,b,ab )
%Law of Cosines
c = sqrt( a^2 + b^2 - 2 * a * b * cosd( ab ));


end

