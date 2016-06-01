function results = getRoughnessVariables(results, sys)

% get the present number of variables
nvars = length(results.iteration.variable);
% add 2 new vars
results.iteration.variable(nvars+1).name = 'Top of biofilm [um]';
results.iteration.variable(nvars+2).name = 'Roughness [um]';
for j = 1:length(results.iteration.variable(1).value),
    iteration = num2str(j);
    fn = [results.directory '\biofilmFront\iteration' iteration '.txt'];
    biomass = load(fn, 'ascii'); 
    [m,n] = size(biomass);
    dx = sys.X / n;
    x = linspace(dx/2, sys.X - dx/2, n);
    y = linspace(sys.X - dx/2, dx/2, m);
    [x,y] = meshgrid(x,y);
    labelBiomass = biomass .* y;
    front = max(labelBiomass, [], 1);
    results.iteration.variable(nvars+1).value(j) = max(front);
    results.iteration.variable(nvars+2).value(j) = std(front);
end;
