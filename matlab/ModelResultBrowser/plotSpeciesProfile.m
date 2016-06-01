function [] = plotSpeciesProfile(results, haxes)

axes(haxes);

% get the system size (TODO later)
sys.X = 500;
sys.Y = sys.X;

% get the concnetration matrix from the file
data = readCurrentSpeciesData(results);

% compute average profile
dataProfile = mean(data,2);
y = linspace(sys.X, 0, length(dataProfile));

hpl = plot(dataProfile, y, 'r-');
set(hpl, 'LineWidth', 1);
hxl = xlabel(['[kg/m^3]']);
set(hxl, 'FontSize', 8);
hyl = ylabel(['distance to surface [\mum]']);
set(hyl, 'FontSize', 8);
set(haxes, 'XLim', [0 max(dataProfile)]);
set(haxes, 'YLim', [0 sys.Y]);
set(haxes, 'FontSize', 5);