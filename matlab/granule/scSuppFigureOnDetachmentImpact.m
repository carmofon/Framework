% scSuppFigureOnDetachmentImpact - to answer reviewers requests

simName(1).name = 'FullAerobicDO100Precision0_01MaxF0_95_eta05';
simName(2).name = 'FullAerobicDO40Precision0_01MaxF0_95_eta05';
simName(3).name = 'AnaerobicFeedDO40Precision0_01MaxF0_95_eta05';
simName(4).name = 'AnaerobicFeedDO20Precision0_01MaxF0_95_eta05';
simName(5).name = 'AnaerobicFeedControlledDO40Precision0_01MaxF0_95_eta05';

figure(1);
set(gcf, 'Position', [428    32   410   890]);
for i = 1:length(simName)
%for i = 1
    subplot(length(simName), 1, i); 
    h = suppFigureOnDetachmentImpactPane(simName(i).name);
    title(sprintf('Simulation #%d', i));
    drawnow;
end;
