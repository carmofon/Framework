% scMakeMovie- Script to make movie frames

dirName = 'D:\results\EPS simulations\TwoSpeciesEpsOxygen';
outDirName = 'D:\lixo\';


results = getResultsFromDirectory(dirName);

results.solutes.current = length(results.solutes.available);

% set the results so that EPS concentration is viewed
results.solids.show = 0;
results.solutes.current = 1;

for i = 1:length(results.iteration.available);
    results.iteration.current = i;
    figure(1);
    clf;
    drawResultsMovie (gca, results);
    title(sprintf('%0.1f h', results.iteration.time(i)));
    drawnow;
    print('-djpeg100', sprintf('%sbiofilm%d.jpg', outDirName, i));
end;