% compute the shear stress to confirm the unit conversion

baseDir = 'C:\Documents and Settings\xavier\Desktop\results\tube';


% construct the results data structure
results = getResultsFromDirectory(baseDir);

viscosity = 1.0e-3 * 1e3* 1e9 * 3600;
tubeRadius = 13000;
gridSide = 65;
sys.X = tubeRadius / (0.5 - 1 / gridSide) * 1.01;
%sys.X = 2.0841e+03
sys.Y = sys.X;
% voxel side
voxelSide = sys.X/gridSide;
% set the iteration to draw
%results.iteration.current = 2;
results.iteration.current = length(results.iteration.variable(1).value)-2;

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = sprintf('%08d', itern); %string

% get either the solute or solid concentration
dn = results.directory;

% load velocity
variable = results.flow.available{2};
fn = [dn '/flow/iteration' iteration '/' variable '.txt'];
v = load(fn, 'ascii');
% load shear
variable = results.flow.available{1};
fn = [dn '/flow/iteration' iteration '/' variable '.txt'];
shearLoaded = load(fn, 'ascii');

figure(1);
set(gcf, 'Position', [74 319 1183 470], 'Color', [1 1 1]);

subplot(1,3,1);
imagesc(v);
axis equal;
colorbar;

% compute the gradient
%compute the gradient
[stressX,stressY] = gradient(-v, voxelSide);
stress = viscosity*sqrt(stressX.^2 + stressY.^2);

subplot(1,3,2);
imagesc(stress);
axis equal;
colorbar;

subplot(1,3,3);
imagesc(shearLoaded);
axis equal;
colorbar;
