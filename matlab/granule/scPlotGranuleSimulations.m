% scPlotGranuleSimulations.m - Plot the main state variables from a
% simulation of 2D granule growing

% add main functions directory to path
addpath('../ModelResultBrowser');

resultPath = 'D:\joao/results';

% the simulation results file to read
%dirName = [resultPath '/FullAerobicDO100Precision0_01MaxF0_95'];
%dirName = [resultPath '/FullAerobicDO40Precision0_01MaxF0_95'];
%dirName = [resultPath '/AnaerobicFeedDO40Precision0_01MaxF0_95'];
%dirName = [resultPath '/AnaerobicFeedDO20Precision0_01MaxF0_99'];
%dirName = [resultPath '/AnaerobicFeedControlledDO40Precision0_01MaxF0_95'];

%dirName = [resultPath '/FullAerobicDO100Test'];
%dirName = [resultPath '/FullAerobicDO100Precision0_01MaxF0_95_Corrected'];
%dirName = [resultPath '/FullAerobicDO40Precision0_01MaxF0_95_Corrected'];
%dirName = [resultPath '/AnaerobicFeedDO40Precision0_01MaxF0_95_Corrected'];
dirName = [resultPath '/AnaerobicFeedDO20Precision0_01MaxF0_95_Corrected'];
%dirName = [resultPath '/AnaerobicFeedControlledDO40Precision0_01MaxF0_95_Corrected'];

% read the results file
results = getResultsFromDirectory(dirName);

% list the variables and their indexes
for i = 1:length(results.iteration.variable),
    str = sprintf('%s = results.iteration.variable(%d).value;',...
        results.iteration.variable(i).name, i);
    disp(str);
end;

% recover the variables
iteration = results.iteration.variable(1).value;
time = results.iteration.variable(2).value;
realTime = results.iteration.variable(3).value;
biovolume = results.iteration.variable(4).value;
runLengthX = results.iteration.variable(5).value;
runLengthY = results.iteration.variable(6).value;
runLengthZ = results.iteration.variable(7).value;
detachedBiomass = results.iteration.variable(8).value;
erodedBiomass = results.iteration.variable(9).value;
sloughedBiomass = results.iteration.variable(10).value;
producedBiomass = results.iteration.variable(11).value;
biomass = results.iteration.variable(12).value;
oxygen = results.iteration.variable(13).value;
oxygenRate = results.iteration.variable(14).value;
ammonium = results.iteration.variable(15).value;
ammoniumRate = results.iteration.variable(16).value;
nitrite = results.iteration.variable(17).value;
nitriteRate = results.iteration.variable(18).value;
nitrate = results.iteration.variable(19).value;
nitrateRate = results.iteration.variable(20).value;
substrate = results.iteration.variable(21).value;
substrateRate = results.iteration.variable(22).value;
phosphate = results.iteration.variable(23).value;
phosphateRate = results.iteration.variable(24).value;
activeNH = results.iteration.variable(25).value;
inertNH = results.iteration.variable(26).value;
activeNO = results.iteration.variable(27).value;
inertNO = results.iteration.variable(28).value;
activeH = results.iteration.variable(29).value;
inertH = results.iteration.variable(30).value;
activePAO = results.iteration.variable(31).value;
phbPAO = results.iteration.variable(32).value;
polypPAO = results.iteration.variable(33).value;
glycogenPAO = results.iteration.variable(34).value;
inertPAO = results.iteration.variable(35).value;

totalInert = inertNH + inertNO + inertH + inertPAO;
timeDay = time/24;

iterationToShowIndex = unique(iteration(find(mod(time, 15) == 0)));


% find the indexes for the begining and end of the cycles
beginCycleTimes = 0:3:time(end);
[uniqueTime, iUT, jUT] = unique(time);


%end indexes
estimatedEndIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes(2:end) - 0.001, 'linear') ;

endOfCycleIndex = round(estimatedEndIndexes);

%begining indexes
estimatedBeginIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes + 0.001, 'linear') ;

beginCycleIndex = ceil(estimatedBeginIndexes);
beginCycleIndex = beginCycleIndex(1:end-1);

% create a command to format each of the subplots
plotFormat = ['set(gca,',...
    ' ''Color'', [0 0 0],',...
    ' ''XLim'', rangeToPlot,',...
    ' ''XColor'', [1 1 1],',...
    ' ''YColor'', [1 1 1])'];

substrateColor = [0 0 1];
oxygenColor = [1 1 1];
ammoniumColor = [0 1 0];
nitriteColor = [1 1 0];
nitrateColor = [0 0.5 0];
phosphateColor = [1 0 0];

activePAOColor = phosphateColor;
phbPAOColor = [0 0.5 0];
polypPAOColor = [1 1 1];
glycogenPAOColor = [0 1 1];
totalInertColor = [0.5 0.5 0.5];
activeNHColor = ammoniumColor;
activeNOColor = nitriteColor;
activeHColor = substrateColor;

%% plot the figures
rangeToPlot = [0 time(end)/24];

% draw the granule run length
figure(1);
set(gcf, 'Color', [0 0 0]);
averageRunLength = (runLengthX + runLengthY) * 0.5;
hl = plot(timeDay, runLengthX, 'w:', timeDay, runLengthY, 'w--',...
    timeDay, averageRunLength, 'w-');
set(hl(3), 'LineWidth', 2);
xlabel('time [day]');
ylabel({'granule diameter';'[{\mu}m]'});
eval(plotFormat);

% draw the particulates
figure(2);
set(gcf, 'Color', [0 0 0]);
h = plot(timeDay, activePAO*1e-9, 'm-');
set(h, 'Color', activePAOColor);
hold on;
h = plot(timeDay, phbPAO*1e-9, 'g-');
set(h, 'Color', phbPAOColor);
h = plot(timeDay, polypPAO*1e-9, 'g-');
set(h, 'Color', polypPAOColor);
h = plot(timeDay, glycogenPAO*1e-9, 'g-');
set(h, 'Color', glycogenPAOColor);
h = plot(timeDay, totalInert*1e-9, 'w-');
set(h, 'Color', totalInertColor);
h = plot(timeDay, activeNH*1e-9, 'b-');
set(h, 'Color', activeNHColor);
h = plot(timeDay, activeNO*1e-9, 'r-');
set(h, 'Color', activeNOColor);
h = plot(timeDay, activeH*1e-9, 'r-');
set(h, 'Color', activeHColor);
hold off;

[legh, objh, outh, outm] = legend(...
    'XPAO',...
    'PHB',...
    'PP',...
    'GLY',...
    'I',...
    'XNH',...
    'XNO',...
    'XH'...
    );
legend('boxoff');
% set the legend properties
htext = text('string', 'Particulates', 'Color', [1 1 1],...
    'Position', [536 4 0]);
set(legh, 'Position', [0.66 0.588 0.122 0.306]);
for i = 1:length(objh),
    if(strmatch(get(objh(i), 'Type'), 'text'))
        set(objh(i), 'Color', [1 1 1]);
    end;
end;
%
xlabel('Time [day]');
ylabel({'Biomass'; '[{\mu}g]'});
eval(plotFormat);


%draw concentration at for a few cycles in the end of simulation
figure(3);
set(gcf, 'Color', [0 0 0]);
h = plot(time, substrate);
set(h, 'Color', substrateColor);
hold on;
h = plot(time, oxygen, 'b-');
set(h, 'Color', oxygenColor);
h = plot(time, ammonium, 'b-');
set(h, 'Color', ammoniumColor);
h = plot(time, nitrite, 'r-');
set(h, 'Color', nitriteColor);
h = plot(time, nitrate, 'g-');
set(h, 'Color', nitrateColor);
h = plot(time, phosphate, 'm-');
set(h, 'Color', phosphateColor);
hold off;
xlabel('Time [h]');
ylabel({'Solute concentration'; '[gCOD/L, gN/L or gP/L]'});
eval(plotFormat);
set(gca,'XLim', rangeToPlot*24)
%set(gca,'XLim', [8002    8023]);

%draw concentration at the end of each cycle (output concentrations)
figure(4);
set(gcf, 'Color', [0 0 0]);
h = plot(timeDay(endOfCycleIndex), substrate(endOfCycleIndex));
set(h, 'Color', substrateColor);
hold on;
h = plot(timeDay(endOfCycleIndex), ammonium(endOfCycleIndex), 'b-');
set(h, 'Color', ammoniumColor);
h = plot(timeDay(endOfCycleIndex), nitrite(endOfCycleIndex), 'r-');
set(h, 'Color', nitriteColor);
h = plot(timeDay(endOfCycleIndex), nitrate(endOfCycleIndex), 'g-');
set(h, 'Color', nitrateColor);
h = plot(timeDay(endOfCycleIndex), phosphate(endOfCycleIndex), 'm-');
set(h, 'Color', phosphateColor);
hold off;
[legh, objh, outh, outm] = legend(...
    'S',...
    'NH_4',...
    'NO_2',...
    'NO_3',...
    'PO_4'...
    );
legend('boxoff');
% set the legend properties
set(legh, 'Position', [0.66 0.131 0.115 0.281]);
for i = 1:length(objh),
    if(strmatch(get(objh(i), 'Type'), 'text'))
        set(objh(i), 'Color', [1 1 1]);
    end;
end;
%
xlabel('Time [day]');
ylabel({'Solute concentration'; '[gCOD/L, gN/L or gP/L]'});
eval(plotFormat);
htext = text('string', 'Solutes', 'Color', [1 1 1], 'Position', [536 0.4 0]);


%% analysis of nitrification / desnitrification
ammoniumRate = results.iteration.variable(16).value * 1.5558e-8/3;
nitriteRate = results.iteration.variable(18).value* 1.5558e-8/3;
nitrateRate = results.iteration.variable(20).value* 1.5558e-8/3;

nitrification = - ammoniumRate;
denitrification = - ammoniumRate - nitriteRate - nitrateRate;

figure(5);
set(gcf, 'Color', [0 0 0]);

h = plot(time, ammoniumRate, 'b-');
set(h, 'Color', ammoniumColor);

hold on;

h = plot(time, nitriteRate, 'r-');
set(h, 'Color', nitriteColor);
h = plot(time, nitrateRate, 'g-');
set(h, 'Color', nitrateColor);

hold off;
xlabel('Time [h]');
ylabel({'Solute Rate'; '[gN/L.h]'});
eval(plotFormat);
set(gca,'XLim', rangeToPlot*24)
%set(gca,'XLim', [8002    8023]);

figure(6);
set(gcf, 'Color', [0 0 0]);

h = plot(time, nitrification, 'b-');
set(h, 'Color', ammoniumColor);

hold on;

h = plot(time, denitrification, 'r-');
set(h, 'Color', nitriteColor);

hold off;
xlabel('Time [h]');
ylabel({'Rate of nitrification (green) and denitrification (yellow)';...
    '[gN/L.h]'});
eval(plotFormat);
set(gca,'XLim', rangeToPlot*24);
%set(gca,'XLim', [8002    8023]);

% simultaneous nitrification/denitrification
snd = denitrification;
snd(denitrification > nitrification) = 0; 
and = denitrification;
and(denitrification < nitrification) = 0; 
figure(7);
set(gcf, 'Color', [0 0 0]);
h = plot(time, snd, 'g-',time, and, 'y-');
xlabel('Time [h]');
ylabel({'SND (green) and AND (yellow)';...
    '[gN/L.h]'});
eval(plotFormat);
%set(gca,'XLim', [8002    8023]);

nRemovalBySnd = sum(snd(1:end-1).*diff(time));
nRemovalByAnd = sum(and(1:end-1).*diff(time));
nRemovalTotal = nRemovalBySnd + nRemovalByAnd;
fSnd = nRemovalBySnd/nRemovalTotal;
fAnd = nRemovalByAnd/nRemovalTotal;
disp(sprintf('N-removal: %0.2f SND, %0.2f AND', fSnd, fAnd)); 



%% Analysis of steady state

%----------------Define the maturatio stage start
startOfMaturation = 240;

%------- Average granule composition
matureStateIndexes = find(time > startOfMaturation*24);
activeNHAverage = mean(activeNH(matureStateIndexes));
activeNOAverage = mean(activeNO(matureStateIndexes));
activeHAverage = mean(activeH(matureStateIndexes));
activePAOAverage = mean(activePAO(matureStateIndexes));
phbPAOAverage = mean(phbPAO(matureStateIndexes));
polypPAOAverage = mean(polypPAO(matureStateIndexes));
glycogenPAOAverage = mean(glycogenPAO(matureStateIndexes));
totalInertAverage = mean(totalInert(matureStateIndexes));

total = activeNHAverage + activeNOAverage + activeHAverage +...
    activePAOAverage + phbPAOAverage + polypPAOAverage + totalInertAverage;

activeNHFraction = activeNHAverage/total;
activeNOFraction = activeNOAverage/total;
activeHFraction = activeHAverage/total;
activePAOFraction = activePAOAverage/total;
phbPAOFraction = phbPAOAverage/total;
polypPAOFraction = polypPAOAverage/total;
glycogenPAOFraction = glycogenPAOAverage/total;
totalInertFraction = totalInertAverage/total;

disp(sprintf('activeNH = %0.1f', activeNHFraction *100));
disp(sprintf('activeNO = %0.1f', activeNOFraction *100));
disp(sprintf('activeH = %0.1f', activeHFraction *100));
disp(sprintf('activePAO = %0.1f', activePAOFraction *100));
disp(sprintf('phbPAO = %0.1f', phbPAOFraction *100));
disp(sprintf('polyPPAO = %0.1f', polypPAOFraction *100));
disp(sprintf('glycogenPPAO = %0.1f', glycogenPAOFraction *100));
disp(sprintf('totalInert = %0.1f', totalInertFraction *100));

%------- Analysis of removal
endOfCycleIndexMature =...
    endOfCycleIndex(endOfCycleIndex > matureStateIndexes(1));
ammoniumEffluent = mean(ammonium(endOfCycleIndexMature));
nitriteEffluent = mean(nitrite(endOfCycleIndexMature));
nitrateEffluent = mean(nitrate(endOfCycleIndexMature));
substrateEffluent = mean(substrate(endOfCycleIndexMature));
phosphateEffluent = mean(phosphate(endOfCycleIndexMature));

ammoniumInfluent = ammonium(1);
substrateInfluent = substrate(1);
phosphateInfluent = phosphate(1);

nRemoval = 1 -(ammoniumEffluent + nitriteEffluent + nitrateEffluent)/...
    ammoniumInfluent;
codRemoval = 1 - (substrateEffluent)/substrateInfluent;
pRemoval = 1 - (phosphateEffluent)/phosphateInfluent;

disp(sprintf('nRemoval = %0.1f', nRemoval *100));
disp(sprintf('++ %0.1f by SND', nRemoval *100*fSnd));
disp(sprintf('++ %0.1f by AND', nRemoval *100*fAnd));
disp(sprintf('codRemoval = %0.1f', codRemoval *100));
disp(sprintf('pRemoval = %0.1f', pRemoval *100));


%%
figure(8);
plot(time/24, (realTime-realTime(1))/1000/3600);
xlabel('time [day]')
ylabel('realTime [h]');

figure(10);
plot(time/24, iteration);
xlabel('time [day]')
ylabel('iteration');
