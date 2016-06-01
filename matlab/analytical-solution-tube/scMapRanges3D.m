% scPlotSolutions


% create the dimensionles groups
N = 50
% c1 = logspace(log10(6.17e-15), log10(0.72), N);
% c2 = logspace(log10(6.67e-2), log10(12.2), N);
% c3 = logspace(log10(2.13e-8), log10(640), N);
c1 = linspace((6.17e-15), (0.72), N);
c2 = linspace((6.67e-2), (12.2), N);
c3 = logspace(log10(2.13e-8), log10(640), N);

[c1, c2, c3] = meshgrid(c1, c2, c3);

r_f = zeros(size(c1));
r_delta = zeros(size(c1));
C_S_surf = zeros(size(c1));
solutionType = zeros(size(c1));

for i = 1:length(c1(:)),
    disp(sprintf('solved %d of %d', i, length(c1(:))));
    % solve
    [r_f(i), r_delta(i), C_S_surf(i), solutionType(i)] =...
        solveSteadyStateDimensionless(c1(i), c2(i), c3(i));
end;

% save workspace
save solutionsLinear;

n=39;
pcolor(c1(:,:,n), c2(:,:,n), solutionType(:,:,n));
%set(gca, 'XScale', 'log', 'YScale', 'log');
colorbar;