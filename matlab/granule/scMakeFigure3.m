% scMakeFigure3.m - Make figure 3 from the paper

% add main functions directory to path
addpath('../ModelResultBrowser');
%resultPath = '/mnt/hda5/people/jxavier/results';
%basedir = '/home/jxavier/';
basedir = 'e:\jxavier';
resultPath = [basedir '/results'];


% the simulation results file to read
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


% defining the colors
substrateColor = [0 0 1];
oxygenColor = [0 0 0];
ammoniumColor = [0 1 0];
nitriteColor = [1 0.8 0];
nitrateColor = [0 0.5 0];
phosphateColor = [1 0 0];

activePAOColor = phosphateColor;
phbPAOColor = [0 0.5 0];
polypPAOColor = [0 0 0];
glycogenPAOColor = [0 1 1];
totalInertColor = [0.5 0.5 0.5];
activeNHColor = ammoniumColor;
activeNOColor = nitriteColor;
activeHColor = substrateColor;

% create a command to format each of the subplots
plotFormat = ['set(gca,''Color'', [1 1 1]),''FontSize'', 8'];
plotAxis = ['set(gca,''XLim'', [0 730])'];



nsims = length(simulation);
figure(1);
set(1, 'Position', [3    10   878   939], 'Color', [1 1 1],...
    'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');

for i = 1:nsims,
    %for i = 1:1,
    %load the data
    results = getResultsFromDirectory(simulation(i).dirName);
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

    %------- draw the particulates
    hparticulate = subplot(nsims, 2, 2*i-1);
    %pos = get(hparticulate, 'Position');
    %set(hparticulate, 'Position', [pos(1) pos(2) pos(3) pos(4)*1.5]);
    %hparticulate = axes('OuterPosition', [0.425 0.441 0.45 0.522]);
    h = plot(timeDay(ia), activePAO(ia)*1e-9, 'm-');
    set(h, 'Color', activePAOColor);
    hold on;
%     h = plot(timeDay(ia), phbPAO(ia)*1e-9, 'g-');
%     set(h, 'Color', phbPAOColor);
%     h = plot(timeDay(ia), polypPAO(ia)*1e-9, 'g-');
%     set(h, 'Color', polypPAOColor);
%     h = plot(timeDay(ia), glycogenPAO(ia)*1e-9, 'g-');
%     set(h, 'Color', glycogenPAOColor);
    h = plot(timeDay(ia), totalInert(ia)*1e-9, 'w-');
    set(h, 'Color', totalInertColor);
    h = plot(timeDay(ia), activeNH(ia)*1e-9, 'b-');
    set(h, 'Color', activeNHColor);
    h = plot(timeDay(ia), activeNO(ia)*1e-9, 'r-');
    set(h, 'Color', activeNOColor);
    h = plot(timeDay(ia), activeH(ia)*1e-9, 'r-');
    set(h, 'Color', activeHColor);
    hold off;
    %
    xlabel('Time [day]');
    set(gca, 'XTick', [0 200 400 600]);
    ylabel({'Particulates';'in granule'; '[{\mu}gCOD]'});
    eval(plotFormat);
    eval(plotAxis);
    set(hparticulate, 'YLim', [0 5]);


    %------- draw concentration at the end of each cycle
    %(output concentrations)
    heffluent = subplot(nsims, 2, 2*i);
    %pos = get(heffluent, 'Position');
    %set(heffluent, 'Position', [pos(1) pos(2) pos(3) pos(4)*1.5]);
    h = plot(timeDay(ia3), substrate(ia3));
    set(h, 'Color', substrateColor);
    hold on;
    h = plot(timeDay(ia3), ammonium(ia3), 'b-');
    set(h, 'Color', ammoniumColor);
    h = plot(timeDay(ia3), nitrite(ia3), 'r-');
    set(h, 'Color', nitriteColor);
    h = plot(timeDay(ia3), nitrate(ia3), 'g-');
    set(h, 'Color', nitrateColor);
    h = plot(timeDay(ia3), phosphate(ia3), 'm-');
    set(h, 'Color', phosphateColor);
    %
    xlabel('Time [day]');
    set(gca, 'XTick', [0 200 400 600]);
    ylabel({'Solute in effluent'; '[gCOD/L,';'gN/L or gP/L]'});
    set(gca, 'YAxisLocation', 'right');
    eval(plotFormat);
    eval(plotAxis);
    set(heffluent, 'YLim', [0 0.1]);
end;

% save the image to file
drawnow;
% %print temporary image with high resolution and save it
% frameName = [basedir '/public_html/transfer/figure3.jpg'];
% frameBig = [basedir '/public_html/transfer/figure3_big.tif'];
% print(1, '-dtiff', '-r300', frameBig);
% % re-open the temporary image
% largeImage = imread(frameBig);
% [mLarge, nLarge, nChannels] = size(largeImage);
% % resize it
% mSmall = 600;
% smallImage = imresize(largeImage, mSmall/mLarge, 'bilinear');
% % save it
% imwrite(smallImage, frameName, 'jpg', 'Quality', 95);
% % remove the temp file
% delete(frameBig);
print('-depsc', ['e:\jxavier\results\temp\Figure3-matlab.eps']);
