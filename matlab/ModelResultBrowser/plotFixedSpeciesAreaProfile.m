function [] = plotFixedSpeciesAreaProfile(results, haxes)

axes(haxes);

% get the system size (TODO later)
sys.X = 500;
sys.Y = sys.X;

% get the concnetration matrix from the file
species = readFixedSpeciesData(results);

% compute average profiles
for i = 1:length(species),
    dataProfile(:,i) = mean(species(i).data,2);
end;

y = (linspace(sys.X, 0, length(dataProfile(:,i))))';
for i = 2:length(species),
    y(:,i) = y(:,1);
end;

% change the order between series
eps = dataProfile(:,1); 
active = dataProfile(:,2); 
inert = dataProfile(:,3); 
dataProfile(:,1) = inert;
dataProfile(:,2) = eps;
dataProfile(:,3) = active;

%draw the area
ha = area(y, dataProfile);
set(ha(1), 'FaceColor', [0.5 0.5 0.5], 'LineStyle', 'none');
set(ha(2), 'FaceColor', 'y', 'LineStyle', 'none');
set(ha(3), 'FaceColor', [0.8 0 0], 'LineStyle', 'none');
view([-90 90]);

hxl = xlabel(['distance to surface [\mum]']);
hyl = ylabel(['concentration [kg/m^3]']);
set(hxl, 'FontSize', 8);
set(hyl, 'FontSize', 8);
set(haxes, 'XLim', [0 sys.X]);
set(haxes, 'YLim', [0 max(dataProfile(:))]);
set(haxes, 'FontSize', 5);
legend('inert', 'eps', 'active mass', 2);

%%%
function species = readFixedSpeciesData(results)

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = num2str(itern);                                             %string

for i = 1:length(results.solids.available),
    variable = results.solids.available{i};
    fn = [results.directory '/solids/iteration' iteration '/' variable '.txt'];
    species(i).data = load(fn, 'ascii');
end;
