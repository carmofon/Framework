% scCase1GranuleEvolution_Figure5

addpath('../ModelResultBrowser');
basedir = 'E:\jxavier';
resultPath = [basedir '/results'];

dirName = [resultPath '/FullAerobicDO100Precision0_01MaxF0_95_CN50'];

results = getResultsFromDirectory(dirName);
% the size of the system side
results.sizeX = 1600;

%
iteration = results.iteration.variable(1).value;
time = results.iteration.variable(2).value;
timeDay = time/24;

%% draw the figure
figure(1);
set(1, 'Position', [4 486 1275 438], 'PaperPositionMode', 'manual');
iterations = logspace(1, 4.7, 8);
for i = 1:length(iterations),
    % set the current iteration
    results.iteration.current = round(iterations(i));
    imageNotFound = true;
    while (imageNotFound)
        %load image
        fn = sprintf('%s/renders/it%08d.png',...
            dirName, iteration(results.iteration.current));
        try
            renderImage = imread(fn);
            imageNotFound = false;
            subplot(2, 4, i);
            image(renderImage);
            axis equal tight;
            %
            t = timeDay(results.iteration.current);
            if (t > 2),
                xlabel(sprintf('t = %d day',...
                    floor(timeDay(results.iteration.current))),...
                    'Color', [0 0 0], 'FontSize', 12);
            else
                xlabel(sprintf('t = %d hour',...
                    round(t*24)),...
                    'Color', [0 0 0], 'FontSize', 12);
            end;
            axis equal, axis tight;
            set(gca, 'XTick', [], 'YTick', []);
            break;
        catch
            results.iteration.current = results.iteration.current + 1;
        end;
    end;
    drawnow;
end;

print('-depsc', [basedir '\results\temp\figure5-matlab.eps']);

