% draw a plot with solute field and balls

%baseDir = 'C:\Documents and Settings\xavier\Desktop\results\tube';
baseDir = 'E:\jxavier\results\TubeTest2';


% construct the results data structure
results = getResultsFromDirectory(baseDir);


%select variable to plot
results.solutes.show = 0;
results.solids.show = 0;
results.flow.show = 1;

scale_velocity = 2.7778e-07; %converts um/h to m/s
max_min_velocity = [0 8];

scale_shear = 1e-12/(3600^2); %converts 1e-15g/um/h to kg/m/s
max_min_shear = [0 0.008];
%max_min_shear = [0 0.025];


%create the subdirectory to save
mkdir(baseDir, 'frames');

%for i = 1:length(results.iteration.variable(1).value)
    %for i = 4456:4456
    %for i = length(results.iteration.variable(1).value):length(results.iteration.variable(1).value)
    for i = 137:137
    try,
        disp(sprintf('creating %d of %d', i, length(results.iteration.variable(1).value)));
        % first change the current iteration
        results.iteration.current = i;
        %get the time
        t = results.iteration.variable(2).value(i);
        figure(1);
        set(gcf, 'Position', [125 319 1029 470], 'Color', [1 1 1],...
            'PaperPositionMode', 'auto');
        %draw velocity
        ha1 = subplot(1,2,1);
        results.flow.current = 2;
        drawResults (gca, results, max_min_velocity, scale_velocity);
        title('Velocity field [mm.s^{-1}]');
        xlabel(sprintf('time %0.1f day', t/24));
        colormap jet;
        %draw shear
        ha2 = subplot(1,2,2);
        results.flow.current = 1;
        %drawResults (gca, results, [], scale_shear);
        hColorBar = drawResults (gca, results, max_min_shear, scale_shear);
        title('Shear stress [N.m^{-2}]');
        xlabel(sprintf('time %0.1f day', t/24));
        set(hColorBar, 'YTick', 0:max_min_shear(end)/4:max_min_shear(end));
        set(ha2, 'ActivePositionProperty', 'position');
        set(ha2, 'Position', [0.57    0.1100    0.2980    0.8150]);
        load shearColorMap;
        colormap(shearColorMap);
        drawnow;
%         %save
%         print('-dpng', '-r300',...
%             sprintf('%s\\frames\\flow%06d.png',baseDir, i));
%         imFile = sprintf('%s\\frames\\flow%06d.tif',baseDir, i);
%         print('-dtiff', '-r300',...
%             imFile);
%         % close the figure
%         % open image and resize it
%         im = imread(imFile);
%         imSmall = imresize(im, 0.3, 'bilinear');
%         imwrite(imSmall,...
%             sprintf('%s\\frames\\flow%06d.png',baseDir, i), 'png');
%         delete(imFile);
%         %
%         close 1;
    catch
        disp(lasterr);
    end
end;