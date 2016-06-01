% scShowLevelSet - show detachment level set with particles

dirName = 'D:\results\detachment\constantconcentration\dfreps\lixo\';

% print figure for each iteration with 

mkdir(dirName, 'levelSet');
outDirName = [dirName '\levelSet\'];
results = getResultsFromDirectory(dirName);


iteration = 400;
results.iteration.current = iteration;
time = results.iteration.variable(2).value(iteration)/24;
figure(1);
%him = axes('position',[.25  .4  .5  .5]);
him = gca;
results.solutes.show = 0;
results.particles = 0;
drawResults (him, results);
xlabel(sprintf('Oxygen concentration at %0.1f day [kg/m^3]', time));
%print
print('-djpeg100', sprintf('%soxygen%04d.jpg', outDirName, iteration));
%clf;
