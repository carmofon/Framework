function [x1mg, x2mg, pdf2d] = drawPdfS(x2, x1, params)

NEST = 100;   %number of estimations
NPTS = 1000;   %number of estimations

xmax = params.mu / params.b;
xEstimate = linspace(0, xmax*1.1, NEST);

pdf2d = zeros(NEST);

for i = 2:length(xEstimate),
    xi = xEstimate(i);
    f = rand(1, NPTS);

    dist1Dt = sdeNextStep(ones(1, NPTS)*xi, params);

    pdf = ksdensity(dist1Dt, xEstimate, 'npoints', 1000, 'kernel', 'normal');

    figure(100);
    plot3(xi*ones(size(xEstimate)), xEstimate, pdf, 'b-');
    hold on;
    xlabel('X_n');
    ylabel('X_{n+1}');
    drawnow;

    pdf2d(:, i) = pdf;
end;

clf(100);

[x1mg, x2mg] = meshgrid(xEstimate, xEstimate);

log10pdf2d = log10(pdf2d);
log10pdf2d(log10pdf2d < -1) = NaN;

figure(2);
pcolor(x1mg, x2mg, log10pdf2d);
colorbar;
shading flat;
