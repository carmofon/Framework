% sc_EstimationForFigure2

% Model
params.mu = 0.5;
params.xMax = 3.5;
params.b = params.mu / params.xMax;
params.alpha = 1000;
params.integrationStepsNumber = 100;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;
params.kernelWidth = 0.1;
x0 = 0.01;

% read the time curve data:
baseDir = 'D:\joao\sloughing-paper';
dirName = [baseDir '\data'];
dirOut = [baseDir '\figs_temp'];
%fileName = '/r01l1i11days.txt';
fileName = '/r03l1i21days.txt';
%
growthCurveData = readGrowthCurve([dirName fileName]);
%
pairs = [];
for i = 1:length(growthCurveData.sensor),
    x1 = growthCurveData.sensor(i).value(1:end-1);
    x2 = growthCurveData.sensor(i).value(2:end);
    if isempty(pairs),
        pairs = [x1'; x2'];
    else
        pairs = [pairs(1,:), x1'; pairs(2,:), x2'];
    end;
end;


% the time series
t = 0:params.samplingStepSize:36;

% extract the pair data
x1 = pairs(1, :);
x2 = pairs(2, :);

figure(4);
hold on;
h = plot(x1, x2, 'b-');
hold off;


% measure the likelihood for the set of parameters
muScanningInterval = logspace(log10(0.1), log10(0.5), 20);
alphaScanningInterval = logspace(log10(100), log10(10000), 20);

[muGrid, alphaGrid] = meshgrid(muScanningInterval, alphaScanningInterval);
L = zeros(size(muGrid));
L(:) = NaN;

i = 0;
for i = 1:length(muScanningInterval);
    for j = 1:length(alphaScanningInterval);
        params.mu = muScanningInterval(i);
        params.xMax = 3.5;
        params.b = params.mu / params.xMax;
        params.alpha = alphaScanningInterval(j);
        [L(j, i), params] = likelihoodOfModelParameters(x2(:), x1(:), params);
        kernelWidth(i) = params.kernelWidth;
        figure(8);
        pcolor(muGrid, alphaGrid, L);
        xlabel('mu');
        ylabel('alpha');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        colorbar;
        drawnow;
    end;
end;

figure(9);
contourf(muGrid, alphaGrid, L, 50);
xlabel('mu');
ylabel('alpha');
set(gca, 'XScale', 'log', 'YScale', 'log');
colorbar;
