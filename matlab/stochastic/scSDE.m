% scSDE

NCURVES = 1000;

x1 = 0.1;

params.mu = 0.5;
params.b = 1/6;
params.alpha = 500;
params.integrationStepsNumber = 100;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;



% solve logistic analtically
tan = 0:params.samplingStepSize/10:36;
xanalitical = params.mu * x1 ./ (params.b * x1 +...
    (params.mu - params.b*x1) * exp (-params.mu * tan));


% sde numerical solutions
t = 0:params.samplingStepSize:36;
x2 = ones(1, NCURVES)*x1;
for i = 2:length(t);
    x2(i,:) = sdeNextStep(x2(i-1, :), params);
end;

meanSde = mean(x2, 2);

figure(1);
h = plot(tan, xanalitical, 'k-', t, x2, 'b-', t, meanSde, 'r-');
set(h(end), 'LineWidth', 2);
xlabel('t');
ylabel('X');
set(gca, 'YLim', [0 3.2]);


