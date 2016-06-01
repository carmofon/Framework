function [] = makeMovie(dirName)

mkdir(dirName, 'movie');
outDirName = [dirName '\movie\'];
results = getResultsFromDirectory(dirName);


for iteration = 1:length(results.iteration.variable(1).value),
    createFrame(results, iteration, outDirName)
end;