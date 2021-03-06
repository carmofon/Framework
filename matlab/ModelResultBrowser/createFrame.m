function [] = createFrame(results, iteration, outDirName)

    results.iteration.current = iteration;
    time = results.iteration.variable(2).value(iteration)/24;
    figure(1);
    clf;
    him = axes('position',[.18  .1  .9  .8]);
    %results.solutes.current = 2;
    drawResults (him, results);
    %xlabel(sprintf('Auto-inducer concentration at %0.1f day [kg/m^3]', time));
    xlabel(sprintf('Oxygen concentration at %0.1f day [kg/m^3]', time));
    % draw other solute on top
%     results.solutes.current = 1;
%     drawConcentrationCurves (him, results);
    
%     % draw oxygen concentration quantity plot
    hpl = axes('position',[.3  .78  .5  .18]);
   % plotVariable(results, 3, 1, hpl);    
    plotVariable(results, 7, -1, hpl);    
    % draw biofilm quantity plot
    %TODO: uncomment to draw biofilm thickness plot
%     hpl = axes('position',[.1  .5  .5  .18]);
%     plotVariable(results, 3, 1, hpl);    
    %print
    hpl = axes('position',[.07  .11  .15  .79]);
    % draw the fixed species profile
    plotFixedSpeciesAreaProfile(results, hpl)
    %plotSpeciesProfile(results, hpl);
    drawnow;
    print('-djpeg100', sprintf('%s/ai%04d.jpg', outDirName, iteration));
