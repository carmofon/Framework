% scRemovalAnalysis.m - Plot the main state variables from a
% simulation of 2D granule growing

% add main functions directory to path
basedir = 'e:\jxavier';
resultPath = [basedir '/results'];


% the simulation results file to read
simulation(1).dirName =...
    [resultPath '/FullAerobicDO100Precision0_01MaxF0_95_CN50'];
simulation(1).matureStage = 150;
simulation(2).dirName =...
    [resultPath '/FullAerobicDO40Precision0_01MaxF0_95_CN50'];
simulation(2).matureStage = 110;
simulation(3).dirName =...
    [resultPath '/AnaerobicFeedDO40Precision0_01MaxF0_95_CN50'];
simulation(3).matureStage = 240;
simulation(4).dirName =...
    [resultPath '/AnaerobicFeedDO20Precision0_01MaxF0_95_CN50'];
simulation(4).matureStage = 240;
simulation(5).dirName =...
    [resultPath '/AnaerobicFeedControlledDO40Precision0_01MaxF0_95_CN50'];
simulation(5).matureStage = 240;

for i = 1:length(simulation),
    [nRemoval, fSnd, fAnd] = nRemovalAnalysis(simulation(i).dirName,...
        simulation(i).matureStage);
    simulation(i).nRemoval = nRemoval;
    simulation(i).fSnd = fSnd;
    simulation(i).fAnd = fAnd;
end;

%% draw the bar plot
nRemovalBarPlot = [];
for i = 1:length(simulation),
    nRemoval = simulation(i).nRemoval;
    fSnd = simulation(i).fSnd;
    fAnd = simulation(i).fAnd;
    nRemovalBarPlot = [nRemovalBarPlot;...
        nRemoval*fSnd*100 nRemoval*fAnd*100];
end;
figure(2);
set(2, 'Position', [360   547   605   377], 'Color', [1 1 1],...
    'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
h = bar(nRemovalBarPlot,'stack');
set(gca, 'YLim', [0, 60]);
% set(gca, 'XTickLabel',...
%      {'1'; '2'; '3'; '4'; '5'});
ylabel({'Average mature-stage' 'N-Removal [%]'});
xlabel('Simulation #');

barWidth = 0.70;
set(h(1), 'FaceColor', [0.3 0.3 0.3], 'BarWidth', barWidth); 
set(h(2), 'FaceColor', 'none', 'BarWidth', barWidth); 
legend('Simultaneous nitrification/denitrification (SND)',...
    'Alternated nitrification/denitrification (AND)', 'Location', 'NorthOutside')

print('-depsc', ['D:\joao\results\temp\Figure6-matlab.eps']);