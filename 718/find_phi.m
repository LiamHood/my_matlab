clear; close all; clc;
syms eps phi 
alpha = 5
eq1 = eps == 90-phi-alpha
eq2 = phi == -eps + acosd((6378.1/(6378.1+579.6350))*cosd(eps))
vpasolve(eq1,eq2)
