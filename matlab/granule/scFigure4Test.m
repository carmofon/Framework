simName = 'FullAerobicDO100Precision0_01MaxF0_95';

% add main functions directory to path
addpath('../ModelResultBrowser');
%resultPath = '/mnt/hda5/people/jxavier/results';
basedir = 'D:\joao';
%basedir = '~/';
resultPath = [basedir '/results'];


% the simulation results file to read
simulation(1).dirName = [resultPath '/' simName];

% defining the colors
substrateColor = [0 0 1];
oxygenColor = [0 0 0];
ammoniumColor = [0 1 0];
nitriteColor = [1 0.8 0];
nitrateColor = [0 0.5 0];
phosphateColor = [1 0 0];

activePAOColor = phosphateColor;
phbPAOColor = [0 0 1];
polypPAOColor = [1 0 0];
glycogenPAOColor = [0 1 1];
totalInertColor = [0.5 0.5 0.5];
activeNHColor = ammoniumColor;
activeNOColor = nitriteColor;
activeHColor = substrateColor;

% create a command to format each of the subplots
plotFormat = ['set(gca,''Color'', [1 1 1]),''FontSize'', 8'];
plotAxis = ['set(gca,''XLim'', [0 730])'];



nsims = length(simulation);
i=1;
%for i = 1:nsims,
%for i = 1:1,
%load the data
results = getResultsFromDirectory(simulation(i).dirName);
% the size of the system side
results.sizeX = 1600;

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
averageRunLength = (runLengthX + runLengthY) * 0.5;


iterationToShowIndex = unique(find(mod(time, 15) == 0));

% find the indexes for the begining and end of the cycles
beginCycleTimes = 0:3:time(end);
[uniqueTime, iUT, jUT] = unique(time);

% cycle indexes:
%end indexes
estimatedEndIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes(2:end) - 0.001, 'linear') ;

endOfCycleIndex = round(estimatedEndIndexes);

%begining indexes
estimatedBeginIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes + 0.001, 'linear') ;

beginCycleIndex = ceil(estimatedBeginIndexes);
beginCycleIndex = beginCycleIndex(1:end-1);

%draw the plots
iterationIndex = length(iterationToShowIndex)-1;
ia = 1:iterationToShowIndex(iterationIndex); %index array
%find the time of the index of the end of the present
%cycle
interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
ia2 = ia(end):interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
% time in end of cycle array
ia3 =...
    endOfCycleIndex(find(time(endOfCycleIndex)<=time(ia(end))+3));

%------- Average granule composition
matureStateIndexes = find(time > 150*24);
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

totalN = nitrate + nitrite + ammonium;

% plot compositions during cycle
beginCycleIndexMature =...
    beginCycleIndex(beginCycleIndex > matureStateIndexes(1));

%% Compile data from mature cycles
timeMatureCycles = [];
substrateMatureCycles = [];
phosphateMatureCycles = [];
ammoniumMatureCycles = [];
nitriteMatureCycles = [];
nitrateMatureCycles = [];
oxygenMatureCycles = [];
phbPAOMatureCycles = [];
polypPAOMatureCycles = [];
glycogenPAOMatureCycles = [];
%for i = 400:400,
%for i = 1:100:(length(beginCycleIndexMature)-1),
for i = (length(beginCycleIndexMature)-1); 
    a = beginCycleIndexMature(i);
    b = endOfCycleIndexMature(i+1);
    c = a:b;
    timeCycle = [time(c) - time(a)];
    timeMatureCycles = [timeMatureCycles timeCycle'];
    substrateMatureCycles = [substrateMatureCycles substrate(c)'];
    phosphateMatureCycles = [phosphateMatureCycles phosphate(c)'];
    ammoniumMatureCycles = [ammoniumMatureCycles ammonium(c)'];
    nitriteMatureCycles = [nitriteMatureCycles nitrite(c)'];
    nitrateMatureCycles = [nitrateMatureCycles nitrate(c)'];
    oxygenMatureCycles = [oxygenMatureCycles oxygen(c)'];
    % and storage compounds
    phbPAOMatureCycles = [phbPAOMatureCycles phbPAO(c)'];
    polypPAOMatureCycles = [polypPAOMatureCycles polypPAO(c)'];
    glycogenPAOMatureCycles = [glycogenPAOMatureCycles glycogenPAO(c)'];
end;

labelTextSize = 14;
textSize = 12;
font = 'Arial';
marker = '-';
lineWidth = 1;

set(gcf, 'Position', [6   186   920   586], 'PaperPositionMode', 'auto');
axes;
h = plot(timeMatureCycles, substrateMatureCycles, marker);
set(h, 'Color', substrateColor, 'LineWidth', lineWidth);
hold on;
h = plot(timeMatureCycles, phosphateMatureCycles, marker);
set(h, 'Color', phosphateColor, 'LineWidth', lineWidth);
h = plot(timeMatureCycles, nitrateMatureCycles, marker);
set(h, 'Color', nitrateColor, 'LineWidth', lineWidth);
h = plot(timeMatureCycles, ammoniumMatureCycles, marker);
set(h, 'Color', ammoniumColor, 'LineWidth', lineWidth);
h = plot(timeMatureCycles, nitriteMatureCycles, marker);
set(h, 'Color', nitriteColor, 'LineWidth', lineWidth);
h = plot(timeMatureCycles, oxygenMatureCycles, marker);
set(h, 'Color', oxygenColor, 'LineWidth', lineWidth);
hold off;
xlabel('Time in cycle [h]', 'FontSize', labelTextSize, 'FontName', font);
ylabel({'Solutes in bulk' '[gCOD.L^{-1}, gN.L^{-1}, gP.L^{-1}]'},...
    'FontSize', labelTextSize, 'FontName', font);
set(gca, 'FontSize', textSize, 'FontName', font);
set(gca, 'Units', 'pixels',...
    'Position', [118.0000  363.0000  405.0000  178.0000]);


%% draw the internal storage compounds
phbPAOMatureCycles2 = phbPAOMatureCycles - phbPAOMatureCycles(1);
polypPAOMatureCycles2 = polypPAOMatureCycles - polypPAOMatureCycles(1);
glycogenPAOMatureCycles2 = glycogenPAOMatureCycles - glycogenPAOMatureCycles(1);
% convert to grams per reactor
phbPAOMatureCycles2 = phbPAOMatureCycles2 * 1e-15;
polypPAOMatureCycles2 = polypPAOMatureCycles2 * 1e-15;
glycogenPAOMatureCycles2 = glycogenPAOMatureCycles2 * 1e-15;



axes;
set(gca, 'Units', 'pixels');
set(gca, 'Position', [118.4645   74.3333  405.5355  178.6667]);
h = plot(timeMatureCycles, phbPAOMatureCycles2, '-');
set(h, 'Color', phbPAOColor, 'LineWidth', lineWidth);
hold on;
h = plot(timeMatureCycles, polypPAOMatureCycles2, '-');
set(h, 'Color', polypPAOColor, 'LineWidth', lineWidth);
h = plot(timeMatureCycles, glycogenPAOMatureCycles2, '-');
set(h, 'Color', glycogenPAOColor, 'LineWidth', lineWidth);
hold off;
set(gca, 'FontSize', textSize, 'FontName', font);
xlabel('Time in cycle [h]', 'FontSize', labelTextSize, 'FontName', font);
ylabel({'Changes in storage' 'compounds [gCOD, gP]'},...
    'FontSize', labelTextSize, 'FontName', font);



%% get the iteration of begining of cycle to plot
cyclesStored = beginCycleIndex(1:400:end);
cycleStart = cyclesStored(end);
cycleEnd = endOfCycleIndex(find(beginCycleIndex == cycleStart));
cycleBeginingIteration = iteration(cycleStart);
timeCycle = time(cycleStart);
timeDay = round(timeCycle/24);
fileIndex = iteration(cycleStart) - 1;

% The granule radial profiles of SOLIDS
results.iteration.current = iteration(cycleStart) + 4;
results.solids.show = 1;
results.solutes.show = 0; 

results.solids.current = 2;
[r, vXNH]= radialProfile(results);
results.solids.current = 3;
[r, vXNO]= radialProfile(results);
results.solids.current = 1;
[r, vXH]= radialProfile(results);
results.solids.current = 4;
[r, vXPAO]= radialProfile(results);
% inerts
vI = [];
for i = 6:9,
    results.solids.current = i;
    [r, tempInert]= radialProfile(results);
    if (isempty(vI))
        vI = tempInert;
    else
        vI = vI + tempInert;
    end;
end;
axes;
h = plot(r, vXNH, '-');
set(h, 'Color', activeNHColor, 'LineWidth', lineWidth);
hold on;
h = plot(r, vXNO, '-');
set(h, 'Color', activeNOColor, 'LineWidth', lineWidth);
h = plot(r, vXH, '-');
set(h, 'Color', activeHColor, 'LineWidth', lineWidth);
h = plot(r, vXPAO, '-');
set(h, 'Color', nitrateColor, 'LineWidth', lineWidth);
h = plot(r, vXPAO, '-');
set(h, 'Color', activePAOColor, 'LineWidth', lineWidth);
h = plot(r, vI, '-');
set(h, 'Color', totalInertColor, 'LineWidth', lineWidth);
hold off;
xlabel({'Distance from granule center' '[{\mu}m]'},...
    'FontSize', labelTextSize, 'FontName', font);
ylabel('Particulates [gCOD.L^{-1}]', 'FontSize', labelTextSize...
    , 'FontName', font);
set(gca, 'FontSize', textSize, 'FontName', font);
set(gca, 'Units', 'pixels');
set(gca, 'Position', [607.0000   78.0000  252.0000  178.0000]);
set(gca, 'XLim', [0 700]);

%% load image
axes;
fn = sprintf('%s/renders/it%08d.png', simulation(1).dirName, fileIndex);
renderImage = imread(fn);
image(renderImage([1:472],[90:570],:));
set(gca, 'XTick', [], 'YTick', []);
axis equal tight;
set(gca, 'Units', 'pixels');
set(gca, 'Position', [590.0000  283.0000  274.0000  259.0000]);

print('-depsc', ['D:\joao\results\temp\' simName '.eps']);

