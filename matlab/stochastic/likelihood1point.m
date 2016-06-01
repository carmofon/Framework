function l = likelihood1point(x2, x1, params)

% determine the conditional density value for x2 given x1 based
% on kernel smoothing of a NPTS number of values obtained from monte-carlo
% simulations

NPTS = 10000;

% monte-carlo the next time step
dist1Dt = sdeNextStep(ones(1, NPTS)*x1, params);

% get likelihood from density estimation using kernel smoothing
l = ksdensity(dist1Dt, x2, 'npoints', 1000, 'kernel', 'normal');