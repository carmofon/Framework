function m = saturationFromFraction(c1, c2, k)

% m = saturationFromFraction(c1, c2, k) 

f = c1./c2;
m =  f ./ (f + k);