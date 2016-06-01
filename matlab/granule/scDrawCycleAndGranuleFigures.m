

simulation(1).name = 'FullAerobicDO100Precision0_01MaxF0_95_CN50';
simulation(2).name = 'FullAerobicDO40Precision0_01MaxF0_95_CN50';
simulation(3).name = 'AnaerobicFeedDO40Precision0_01MaxF0_95_CN50';
simulation(4).name = 'AnaerobicFeedDO20Precision0_01MaxF0_95_CN50';
simulation(5).name = 'AnaerobicFeedControlledDO40Precision0_01MaxF0_95_CN50';

for i = 1:length(simulation),
%for i = 1,
    figure(i);
    createGranuleAndCycleFigure(simulation(i).name);
end;