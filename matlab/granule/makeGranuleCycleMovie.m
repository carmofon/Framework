function [] = makeGranuleCycleMovie(dirName)

%make the directory to store the frames
mkdir(dirName, 'frames');

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
averageRunLength = (runLengthX + runLengthY) * 0.5;


iterationToShowIndex = unique(find(mod(time, 15) == 0));

% create a command to format each of the subplots
plotFormat = ['set(gca,',...
    ' ''Color'', [0 0 0],',...
    ' ''XColor'', [1 1 1],',...
    ' ''YColor'', [1 1 1])'];
plotAxis = ['set(gca,',...
    ' ''XLim'', [0 730])'];

substrateColor = [0 0 1];
oxygenColor = [1 1 1];
ammoniumColor = [0 1 0];
nitriteColor = [1 1 0];
nitrateColor = [0.5 0.5 0.5];
phosphateColor = [1 0 0];

activePAOColor = phosphateColor;
phbPAOColor = [0 0.5 0];
polypPAOColor = [1 1 1];
glycogenPAOColor = [0 1 1];
totalInertColor = [0.5 0.5 0.5];
activeNHColor = ammoniumColor;
activeNOColor = nitriteColor;
activeHColor = substrateColor;


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


for iterationIndex = 1:(length(iterationToShowIndex)-1),
    %for iterationIndex = [1, length(iterationToShowIndex)],
    %for iterationIndex = [length(iterationToShowIndex)-1]
    % errorCorrection(1) = find(iterationToShowIndex==10211);
    % errorCorrection(2) = find(iterationToShowIndex==10268);
    % errorCorrection(3) = find(iterationToShowIndex==10328);
    % errorCorrection(4) = find(iterationToShowIndex==10390);
    % for iterationIndex = errorCorrection
    try
        % try to create the image of the right size
        while (1),
            close all;
            %iterationIndex = iterationToShowIndex(end);
            ia = 1:iterationToShowIndex(iterationIndex); %index array
            %find the time of the index of the end of the present
            %cycle
            interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
            ia2 = ia(end):interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
            % time in end of cycle array
            ia3 =...
                endOfCycleIndex(find(time(endOfCycleIndex)<=time(ia(end))+3));

            %create the figure window
            figure(1);
            set(1, 'Position', [90   375   620   300], 'Color', [0 0 0],...
                'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');

            % draw the render
            harender = axes('Position', [-0.08 -0.1 0.65 1.2]);
            %load image
            fn = sprintf('%s/renders/it%08d.png',...
                dirName, iteration(ia(end)));
            try
                renderImage = imread(fn);
            catch
                renderImage = imread('imageNotFound.bmp');
            end;
            image(renderImage);
            axis equal tight;
            text('string', sprintf('%d days', floor(timeDay(ia(end)))),...
                'Color', [1 1 1], 'Position', [289 13.2 0]);


            % draw the particulates
            hcycle = axes('Position', [0.5597 0.1100 0.4274 0.8667]);
            if and(~isempty(endOfCycleIndex), length(endOfCycleIndex)>1),
                % scale time to time in cycle 0 to 3 h
                timeInCycle = time(ia2) - time(ia2(1));
                h = plot(timeInCycle, substrate(ia2));
                set(h, 'Color', substrateColor, 'LineWidth', 2);
                hold on;
                h = plot(timeInCycle, oxygen(ia2), 'b-');
                set(h, 'Color', oxygenColor, 'LineWidth', 2);
                h = plot(timeInCycle, ammonium(ia2), 'b-');
                set(h, 'Color', ammoniumColor, 'LineWidth', 2);
                h = plot(timeInCycle, nitrite(ia2), 'r-');
                set(h, 'Color', nitriteColor, 'LineWidth', 2);
                h = plot(timeInCycle, nitrate(ia2), 'g-');
                set(h, 'Color', nitrateColor, 'LineWidth', 2);
                h = plot(timeInCycle, phosphate(ia2), 'm-');
                set(h, 'Color', phosphateColor, 'LineWidth', 2);
                hold off;
            end;
            xlabel('Time in cycle [h]');
            ylabel({'Solute concentration'; '[gCOD/L, gN/L or gP/L]'});
            eval(plotFormat);
            %%%LEGEND
            set(hcycle, 'XLim', [0 3]);
            [legh, objh, outh, outm] = legend(...
                'S',...
                'O_2',...
                'NH_4',...
                'NO_2',...
                'NO_3',...
                'PO_4'...
                );
            legend('boxoff','Location','NorthEast');
            % set the legend properties
            %htext = text('string', 'Solutes:', 'Color', [1 1 1],...
            %    'Position', [3.3, 0.2, 0]);
            %set(legh, 'Position', [0.85 0.588 0.122 0.306]);
            for i = 1:length(objh),
                if(strmatch(get(objh(i), 'Type'), 'text'))
                    set(objh(i), 'Color', [1 1 1]);
                end;
            end;
            lixo = 0;
            %
            drawnow;
            %print temporary image with high resolution and save it
            frameName = sprintf('%s/frames/it%08d.jpg',...
                dirName, iterationToShowIndex(iterationIndex));
            frameBig = sprintf('%s/frames/big%08d.tif',...
                dirName, iterationToShowIndex(iterationIndex));
            print(1, '-dtiff', '-r300', frameBig);
            % re-open the temporary image
            largeImage = imread(frameBig);
            [mLarge, nLarge, nChannels] = size(largeImage);
            % resize it
            mSmall = 358;
            smallImage = imresize(largeImage, mSmall/mLarge, 'bilinear');
            % save it
            imwrite(smallImage, frameName, 'jpg', 'Quality', 95);
            % remove the temp file
            delete(frameBig);
            if (iterationIndex == 1),
                infoFirst = imfinfo(frameName);
                %manuanly cheack size
                if and((infoFirst.Width == 740),...
                        (infoFirst.Height == 358)),
                    disp(sprintf('image %d is done', iterationIndex));
                    break;
                else,
                    disp(sprintf('size of image is %d x %d ->try again',...
                        infoFirst.Width, infoFirst.Height));
                end;
            %check if size is different from that of the first image
            else
                infoThis = imfinfo(frameName);
                if and((infoThis.Width == infoFirst.Width),...
                        (infoThis.Height == infoFirst.Height)),
                    % if it is, break, otherwise, try again (draw, print,
                    % resize, and save
                    break;
                else
                    disp(['image ' frameName ' will be redrawn']);
                end;
            end;
        end;
    catch
        warning(sprintf('frame %d not created', iterationIndex));
        disp('The error:');
        disp(lasterr);
    end
end;