function [pdf, xf] = probabilityDesityFucntion(xi, params);

NPTS = 10000;

f = rand(1, NPTS);

dist1Dt = sdeNextStep(ones(1, NPTS)*xi, params);

[pdf, xf] = ksdensity(dist1Dt, 'npoints',...
    1000, 'kernel', 'normal', 'width', params.kernelWidth);