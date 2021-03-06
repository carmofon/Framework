function h = suppFigureOnDetachmentImpactPane(simName)


% add main functions directory to path
addpath('../ModelResultBrowser');

basedir = 'z:/joao/granule-stuff/simulation-results/';

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
results =[];
results.directory = simulation(i).dirName;
results =...
    readSimulationResults([simulation(i).dirName...
    '/simulationResults.txt'], results);
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
ia3 = endOfCycleIndex(find(time(endOfCycleIndex)<=time(ia(end))+3));

%------- Average granule composition
matureStateIndexes  = find(time > 150*24);
activeNHAverage     = mean(activeNH(matureStateIndexes));
activeNOAverage     = mean(activeNO(matureStateIndexes));
activeHAverage      = mean(activeH(matureStateIndexes));
activePAOAverage    = mean(activePAO(matureStateIndexes));
phbPAOAverage       = mean(phbPAO(matureStateIndexes));
polypPAOAverage     = mean(polypPAO(matureStateIndexes));
glycogenPAOAverage  = mean(glycogenPAO(matureStateIndexes));
totalInertAverage   = mean(totalInert(matureStateIndexes));

total = activeNHAverage + activeNOAverage + activeHAverage +...
    activePAOAverage + phbPAOAverage + polypPAOAverage + totalInertAverage;

activeNHFraction    = activeNHAverage/total;
activeNOFraction    = activeNOAverage/total;
activeHFraction     = activeHAverage/total;
activePAOFraction   = activePAOAverage/total;
phbPAOFraction      = phbPAOAverage/total;
polypPAOFraction    = polypPAOAverage/total;
glycogenPAOFraction = glycogenPAOAverage/total;
totalInertFraction  = totalInertAverage/total;

%------- Analysis of removal
endOfCycleIndexMature =...
    endOfCycleIndex(endOfCycleIndex > matureStateIndexes(1));

% plot compositions during cycle
beginCycleIndexMature =...
    beginCycleIndex(beginCycleIndex > matureStateIndexes(1));

%% Compile data from mature cycles
timeMatureCycles = [];
totalMass = [];
totalDetached = [];

totalMassThis = totalInert + activeNH + activeNO + activeH + activePAO +...
    phbPAO + polypPAO + glycogenPAO;

for i = 1:100:(length(beginCycleIndexMature)-1),
%for i = (length(beginCycleIndexMature)-1); 
    a = beginCycleIndexMature(i);
    b = endOfCycleIndexMature(i+1);
    c = a:b;
    timeCycle = [time(c) - time(a)];
    timeMatureCycles = [timeMatureCycles timeCycle'];
    totalMass = [totalMass totalMassThis(c)'];
    totalDetached = [totalDetached cumsum(detachedBiomass(c))'];
end;


%% analyze detached biomass

h = plot(timeMatureCycles, (totalDetached./totalMass)*100, 'k-');
xlabel('time in cycle [h]');
set(gca, 'YLim', [0 1]);
ylabel('% detached');