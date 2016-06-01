function  f = objfun(x, xNext, xPresent, params)

% the function to minimize, i.e. -L

params.mu    = x(1);
params.b     = x(2);
params.alpha = x(3);

f = -likelihoodOfModelParameters(xNext, xPresent, params);