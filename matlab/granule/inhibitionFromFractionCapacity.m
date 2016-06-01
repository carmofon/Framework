function m = inhibitionFromFractionCapacity(c1, c2, k, fMax)

% m = inhibitionFromFractionCapacity(c1, c2, k, fMax) 

f = c1./c2;
m =  (f - fMax) ./ ((f - fMax) + k);