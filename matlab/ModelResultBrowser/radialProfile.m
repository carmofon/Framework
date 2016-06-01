function [r, v, vMax, vMin] = radialProfile(results)

% [v, r] = radialProfile(results) compute the radial profile of variable


%read data from current variable
% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = sprintf('%08d', itern); %string

% get either the solute or solid concentration
fn = results.directory;

if (results.solids.show);
    variable = results.solids.available{results.solids.current};
    fn = [fn '/solids/iteration' iteration '/' variable '.txt'];
else (results.solutes.show);
    variable = results.solutes.available{results.solutes.current};
    fn = [fn '/solutes/iteration' iteration '/' variable '.txt'];
end;

%
data2D = load(fn, 'ascii');

[m, n] = size(data2D);
x = ((0 : n-1)/(n-1) - 0.5)* results.sizeX;
y = x;
[x, y] = meshgrid(x, y);

%
nThetas = 50; % number of steps in the theta direction
nRs = 50; % number of steps in the r direction
r = linspace(9, results.sizeX/2, nRs);
vPrelim = [];
c = 1;
for theta = linspace(0, 2*pi, nThetas),
    xLine = r.*cos(theta);
    yLine = r.*sin(theta);
    vPrelim(c,:) = interp2(x, y, data2D, xLine, yLine, 'linear');
    c = c+1;
end;

v = mean(vPrelim);
vMax = max(vPrelim);
vMin = min(vPrelim);





