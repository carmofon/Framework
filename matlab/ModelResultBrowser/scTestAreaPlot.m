% scTestAreaPlot

basedir = 'D:\results\';

data(1).dirName = [basedir 'drg-detachment-AI-eps\umax5.47E-01kdet1.00E-03rai5.00E-02kai3.50E-04\'];

results = getResultsFromDirectory(data(1).dirName);

figure(1);
plotFixedSpeciesAreaProfile(results, gca);