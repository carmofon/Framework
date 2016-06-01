% scDistributionsOfSDE - plot some relevant distributions of the SDE

NPTS = 10000; % number of points for monte-carlo simulations
NEST = 100;   %number of estimations

ALPH = 200;
x0 = 0.1;

params.mu = 0.5;
params.b = 1/6;
params.alpha = ALPH;
params.integrationStepsNumber = 200;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 2;

xmax = params.mu / params.b;
xEstimate = linspace(0, xmax*1.1, NEST);

pdf2d = zeros(NEST);

for i = 2:length(xEstimate),
    xi = xEstimate(i);
    f = rand(1, NPTS);
    dist1dt = f.*(1 - f).^ALPH;
    dist1Dt = sdeNextStep(ones(1, NPTS)*xi, params);

    pdf = ksdensity(dist1Dt, xEstimate, 'npoints', 1000, 'kernel', 'normal');

    figure(1);
    plot3(xi*ones(size(xEstimate)), xEstimate, pdf, 'b-');
    hold on;
    drawnow;
    xlabel('X_n');
    ylabel('X_{n+1}');

    pdf2d(:, i) = pdf;
end;
 
[x1mg, x2mg] = meshgrid(xEstimate, xEstimate);

log10pdf2d = log10(pdf2d);
log10pdf2d(log10pdf2d < -1) = NaN;

figure(2);
plot3(x1mg, x2mg, log10pdf2d);
xlabel('X_n');
ylabel('X_{n+1}');

figure(3);
plot3(x1mg, x2mg, pdf2d);
xlabel('X_n');
ylabel('X_{n+1}');

%%
figure(4);
pcolor(x1mg, x2mg, log10pdf2d);
colorbar;
shading flat;