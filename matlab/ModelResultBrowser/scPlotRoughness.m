% scPlotRoughness

% get the system size (TODO later)
sys.X = 500;
sys.Y = sys.X;

basedir = 'D:\results\drg\';

data(1).dirName = [basedir 'umax5.47E00'];
data(2).dirName = [basedir 'umax5.47E-01'];
data(3).dirName = [basedir 'umax5.47E-02'];
data(4).dirName = [basedir 'umax5.47E-03'];

for i = 1:length(data),
    results = getResultsFromDirectory(data(i).dirName);    
    results = getRoughnessVariables(results, sys);
    plot(results.iteration.variable(8).value, results.iteration.variable(9).value);
    hold on;
    drawnow;
end;