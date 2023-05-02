syms eps phi 
alpha = 4
eq1 = eps == 90-phi-alpha
eq2 = phi == -eps + acosd((6378.1/(6378.1+554.283473))*cosd(eps))
vpasolve(eq1,eq2)
