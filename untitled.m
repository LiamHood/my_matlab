syms eps phi 
alpha = 25
eq1 = eps == 90-phi-alpha
eq2 = phi == -eps + acosd((6378.1/6678.1)*cosd(eps))
vpasolve(eq1,eq2)
