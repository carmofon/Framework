dirName = 'D:\results\detachment\constantconcentration\dfr\mushroopShapes';
outDirName = [dirName '\movie\'];
results = getResultsFromDirectory(dirName);


for iteration = 1:length(results.iteration.variable(1).value),
    results.iteration.current = iteration;
    time = results.iteration.variable(2).value(iteration)/24;
    figure(1);
    %him = axes('position',[.25  .4  .5  .5]);
    him = gca;
    drawResults (him, results);
    xlabel(sprintf('Oxygen concentration at %0.1f day [kg/m^3]', time));
    % draw oxygen concentration quantity plot
    hpl = axes('position',[.1  .78  .5  .18]);
    plotVariable(results, 7, -1, hpl);    
    drawnow;
    % draw biofilm quantity plot
    hpl = axes('position',[.1  .5  .5  .18]);
    plotVariable(results, 3, 1, hpl);    
    drawnow;
    %print
    print('-djpeg100', sprintf('%soxygen%04d.jpg', outDirName, iteration));
    clf;
end;