d = [1.53 1.57 1.54 1.54 1.50 1.51 1.55 1.54 1.56 1.53];
d = data;
dev = d - mean(d) ;
var_top = 0 ;
for ii = 1:length(d)
    var_top = dev(ii)^2 + var_top ;
end
var = var_top / ( length(d) - 1 ) ;
vara = var_top / length(d) ;
k = length(d)-1 ;
