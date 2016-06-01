% scDrawFig2_SDE

% Model
params.mu = 0.2763;
params.xMax = 3.5;
params.b = params.mu / params.xMax;
params.alpha = 1.1288e+003;
params.integrationStepsNumber = 100;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;
sigma = 0.1;

% analytical solution of model
x1 = 0.01;
t = 0:params.samplingStepSize:36;
xanalitical = params.mu * x1 ./ (params.b * x1 +...
    (params.mu - params.b*x1) * exp (-params.mu * t));

% pure logistic
figure(1);
set(1, 'Position', [40 420 1202 515], 'PaperPositionMode', 'auto');

% phase-space diagram of pure logistic
subplot(1, 2, 2);
h = plot(xanalitical(1:end-1), xanalitical(2:end), 'k--');
set(h, 'LineWidth', 2);

% SDE numeric solution of model
NCURVES = 9;
xSde = ones(1, NCURVES)*x1;
xSdeObserved = ones(1, NCURVES)*x1;
xSdeObserved(1,:) =  xSde(1,:) + randn(1, NCURVES)...
    .* sigma ;
for i = 2:length(t);
    xSde(i,:) = sdeNextStep(xSde(i-1, :), params);
    xSdeObserved(i,:) =  xSde(i,:) + randn(1, NCURVES)...
        .* sigma ;
end;

%truncate xSdeObserved values
xSdeObserved(xSdeObserved < 0) = 0;

for i = 1:NCURVES,
    subplot(3, 6, floor((i-1)/3)*6 + mod(i-1,3) + 1);
    h = plot(t, xanalitical, 'k--', t, xSdeObserved(:,i), 'k-');
    set(h(1), 'LineWidth', 2);
    xlabel('time [day]');
    ylabel('X');
    set(gca, 'XLim', [0, 40], 'YLim', [0, 4]);
end;

subplot(1, 2, 2);
hold on;
h = plot(xSdeObserved(1:end-1, :), xSdeObserved(2:end, :), 'k-');
hold off;
set(h, 'LineWidth', 1);
set(gca, 'XLim', [0, 4], 'YLim', [0, 4]);