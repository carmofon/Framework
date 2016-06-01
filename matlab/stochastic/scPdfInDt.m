% scPdfInDt - the P.D.F. computed analytically
% for an integration time step

NGRD = 100;

% Model
params.mu = 0.9;
params.b = 1/4;
params.alpha = 1.1;
params.integrationStepsNumber = 1;
params.integrationStepSize = 1/params.integrationStepsNumber;
params.samplingStepSize = 1;

x0 = 0.1;
xMax = params.mu / params.b;


xPresent = linspace(x0, xMax*1.1, NGRD);
xNext = linspace(0, xMax*1.1, NGRD);

[xPresent, xNext] = meshgrid(xPresent, xNext);

Pr = 1 - [1 + params.integrationStepSize .*...
    (params.mu - xPresent .* params.b) - xNext ./ xPresent]...
    .^(1./params.alpha);
f = 1/params.alpha./xPresent .* [1 + params.integrationStepSize .*...
    (params.mu - xPresent .* params.b) - xNext ./ xPresent]...
    .^(1./params.alpha - 1);

Pr(imag(Pr) ~= 0) = NaN;
Pr(Pr > 1) = NaN;
f(or(or(imag(f) ~= 0, Pr > 1), imag(Pr) ~= 0)) = NaN;

figure(1);
pcolor(xPresent, xNext, f);
set(1, 'Renderer', 'zbuffer');
xlabel('X_n');
ylabel('X_{n+1}');
title('dX_dXgrowth');
shading interp;
colorbar;

