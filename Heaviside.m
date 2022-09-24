function [H] = Heaviside(phi,alpha,epsilon)
H = 3*(1-alpha)/4*(phi/epsilon-phi.^3/(3*(epsilon)^3)) + (1+alpha)/2;
H(phi>epsilon) = 1;
H(phi<-epsilon) = alpha;
end