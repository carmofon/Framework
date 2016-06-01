% scLikelihoodFields - compute L for a range of mu and alpha

% Model
params.mu = 0.5;
params.xMax = 3.5;
params.b = params.mu / params.xMax;
params.alpha = 1000;
params.integrationStepsNumber = 100;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;
sigma = 0.1;
x0 = 0.01;

% the time series
t = 0:params.samplingStepSize:36;

% SDE numeric solution of model
NCURVES = 9;
xSde = ones(1, NCURVES)*x0;
for i = 2:length(t);
    xSde(i,:) = sdeNextStep(xSde(i-1, :), params);
end;

x1 = xSde(1:end-1, :);
x2 = xSde(2:end, :);

figure(3);
hold on;
h = plot(t, xSde, 'b-');
hold off;
set(h, 'LineWidth', 2);

figure(4);
hold on;
h = plot(x1, x2, 'b-');
hold off;


% measure the likelihood for the set of parameters
%params.mu = 10;
params.kernelWidth = 0.01;
%[L, params] = likelihoodOfModelParameters(x2(:), x1(:), params);

muScanningInterval = logspace(log10(0.2), log10(1), 5);
alphaScanningInterval = logspace(log10(500), log10(2000), 5);

[muGrid, alphaGrid] = meshgrid(muScanningInterval, alphaScanningInterval);

i = 0;
for i = 1:length(muScanningInterval);
    for j = 1:length(alphaScanningInterval);
        clear kernelWidth L;
        close all;
        params.mu = muScanningInterval(i);
        params.xMax = 3.5;
        params.b = params.mu / params.xMax;
        params.alpha = alphaScanningInterval(j);
        for w = logspace(-4, 1, 20),
            params.kernelWidth = w;
            [L(i), params] = likelihoodOfModelParameters(x2(:), x1(:), params);
            kernelWidth(i) = params.kernelWidth;
            figure(7)
            plot(kernelWidth, L, 'b+');
            xlabel('Kernel width');
            ylabel('Likelihood');
            drawnow;
        end;
        indexMax = find(L == max(L));
        kernelWidthGrid(j, i) = kernelWidth(indexMax(1));
        L(j, i) = L(indexMax(1));
    end;
end;

figure(8);
subplot(1, 2, 1);
pcolor(muGrid, alphaGrid, kernelWidthGrid);
xlabel('mu');
xlabel('kernelWidthGrid');
subplot(1, 2, 2);
pcolor(muGrid, alphaGrid, L);
xlabel('mu');
xlabel('kernelWidthGrid');



