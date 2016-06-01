% scDrawFig1_and_Fig3

baseDir = 'D:\joao\sloughing-paper';
dirName = [baseDir '\data'];
dirOut = [baseDir '\figs_temp\'];
%fileName = '/r01l1i11days.txt';
fileName = '/r03l1i21days.txt';

growthCurveData = readGrowthCurve([dirName fileName]);


figure(1);
set(1, 'Position', [40 470 750 465], 'PaperPositionMode', 'auto');
for i = 1:9,
    subplot(3, 3, i);
    plot(growthCurveData.time, growthCurveData.sensor(i).value, 'k-');
    xlabel('time [day]');
    ylabel('X');
    set(gca, 'XLim', [0 40], 'YLim', [0 4]);
end;

print('-depsc', [dirOut 'fig1.eps']); 


% PHASE-SPACE (Experimental data)
% group the data into [X_{N+1}; X_{N}] pairs
figure(3);
set(3, 'Position', [40 655 669 269], 'PaperPositionMode', 'auto');
subplot(1, 2, 2);
pairs = [];
for i = 1:length(growthCurveData.sensor),
    x1 = growthCurveData.sensor(i).value(1:end-1);
    x2 = growthCurveData.sensor(i).value(2:end);
    if isempty(pairs),
        pairs = [x1'; x2'];
    else
        pairs = [pairs(1,:), x1'; pairs(2,:), x2'];
    end;
    %h = plot(x1, x2, 'k-');
    h = plot(x1, x2, 'k+');
    hold on;
end;
hold off;
xlabel('X_n');
ylabel('X_{n+1}');
set(gca, 'XLim', [0 4], 'YLim', [0 4]);


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
figure(3);
subplot(1, 2, 1);
h = plot(t, xanalitical, 'k--');
set(h, 'LineWidth', 2);
set(gca, 'XLim', [0 40], 'YLim', [0 4]);
% simulation
% SDE numeric solution of model
NCURVES = 1;
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
% plot them on the "pure" logistic plot
for i = 1:NCURVES,
    hold on;
    h = plot(t, xSdeObserved(:,i), 'k-');
end;


% phase-space diagram of pure logistic
subplot(1, 2, 2);
hold on;
h = plot(xanalitical(1:end-1), xanalitical(2:end), 'k--');
set(h, 'LineWidth', 2);


print('-depsc', [dirOut 'fig3.eps']); 
