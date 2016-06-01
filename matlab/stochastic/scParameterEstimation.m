% scParameterEstimation - estimate parameters from a phoobia growth curve

dirName = 'E:\jxavier\Sloughing detection\data';
fileName = '/r01l1i11days.txt';

growthCurveData = readGrowthCurve([dirName fileName]);

seriesColors = jet(length(growthCurveData.sensor));

% plot the original time series
figure(1);
title('original data');
for i = 1:length(growthCurveData.sensor),
    h = plot(growthCurveData.time, growthCurveData.sensor(i).value);
    set(h, 'Color', seriesColors(i, :));
    hold on;
end;
hold off;
xlabel('time [day]');
ylabel('value');

% PHASE-SPACE
% group the data into [X_{N+1}; X_{N}] pairs
figure(2);
pairs = [];
for i = 1:length(growthCurveData.sensor),
    x1 = growthCurveData.sensor(i).value(1:end-1);
    x2 = growthCurveData.sensor(i).value(2:end);
    if isempty(pairs),
        pairs = [x1'; x2'];
    else
        pairs = [pairs(1,:), x1'; pairs(2,:), x2'];
    end;
    h = plot(x1, x2, '+');
    set(h, 'Color', seriesColors(i, :));
    hold on;
end;
hold off;
xlabel('X_N');
ylabel('X_{N+1}');

% Model
params.mu = 0.9;
params.b = 1/4;
params.alpha = 400;
params.integrationStepsNumber = 100;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;

% analytical solution of model
x1 = 0.1;
t = 0:params.samplingStepSize:36;
xanalitical = params.mu * x1 ./ (params.b * x1 +...
    (params.mu - params.b*x1) * exp (-params.mu * t));

figure(1);
hold on;
h = plot(t, xanalitical, 'r-');
hold off;
set(h, 'LineWidth', 2);

figure(2);
hold on;
h = plot(xanalitical(1:end-1), xanalitical(2:end), 'r-');
hold off;
set(h, 'LineWidth', 2);

% SDE numeric solution of model
NCURVES = 9;
xSde = ones(1, NCURVES)*x1;
for i = 2:length(t);
    xSde(i,:) = sdeNextStep(xSde(i-1, :), params);
end;

figure(3);
hold on;
h = plot(t, xSde, 'b-');
hold off;
set(h, 'LineWidth', 2);

figure(4);
hold on;
h = plot(xSde(1:end-1, :), xSde(2:end, :), 'b+');
hold off;


x1 = xSde(1:end-1, 1);
x2 = xSde(2:end, 1);


% measure the likelihood for the set of parameters
params.mu = 10;
params.kernelWidth = 0.001;
[L, params] = likelihoodOfModelParameters(x2, x1, params);


%%%%%%%%%%%%%
params.kernelWidth = 0.001;
xNext    = x2;
xPresent = x1;
fToMinimize = @(x) objfun(x, xNext, xPresent, params);

p(1) = params.mu * 1.1;
p(2) = params.b * 0.9;
p(3) = params.alpha * 1.1;

fToMinimize(p)

p0 = p;    % Starting guess
options = optimset('LargeScale','on', 'Display', 'iter');
[p, fval, exitflag, output] = fminunc(fToMinimize, p0, options);
%%%%%%%%%%%%%%