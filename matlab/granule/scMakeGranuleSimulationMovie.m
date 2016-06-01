% scMakeGranuleSimulationMovie.m - Make a movi with the frames
% with plots and 2D granule renders

% add main functions directory to path
addpath('../ModelResultBrowser');
resultPath = '/mnt/hda5/people/jxavier/results';
%resultPath = '~/results';


% the simulation results file to read
movie(1).dirName = [resultPath...
    '/FullAerobicDO100Precision0_01MaxF0_95_Corrected'];
movie(2).dirName = [resultPath...
    '/FullAerobicDO40Precision0_01MaxF0_95_Corrected'];

for i = 1:length(movie)
    makeGranuleCycleMovie(movie(i).dirName);
end;