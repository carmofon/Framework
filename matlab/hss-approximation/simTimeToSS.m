%{
This is a code that calls numericalSim.m repeatedly while varying the
parameters to see if the time to the same fixed-point solution of
f can vary with different parameter multipliers.

Time step length must be specified in addition to start and end times
or won't work.
%}

function fs = simTimeToSS(tspan)

lcXs = [0.5, 0.5; 1, 1; 2, 2];
fs = zeros(length(tspan), 3);

for i=1:size(lcXs, 1)
    lambdaX = lcXs(i, 1);
    cX = lcXs(i, 2);
    [t, y] = numericalSim(tspan, lambdaX, cX);
    fs(:, i) = y(:, 1);
end

figure; hold all

for i=1:size(fs, 2)
    plot(tspan, fs(:, i))
end

legend('lambda = 0.5x, cRs = 0.5x', 'lambda = 1x, cRs = 1x', ...
    'lambda = 2x, cRs = 2x')