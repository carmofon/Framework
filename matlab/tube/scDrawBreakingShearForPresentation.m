% draw a plot with solute field and balls

resDir = '/home/jxavier/work/results';
%resDir = 'E:\jxavier\results\';
res(1).baseDir = [resDir '/TubeTest_ShearInducedDetachment5e16'];
res(2).baseDir = [resDir '/TubeTest_ShearInducedDetachment2e17'];


for d = 1:length(res),
    baseDir = res(d).baseDir;
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

    for i = 1:length(results.iteration.variable(1).value)
    %for i = 539,
        try,
            disp(sprintf('creating %d of %d', i, length(results.iteration.variable(1).value)));
            % first change the current iteration
            results.iteration.current = i;
            %get the time
            t = results.iteration.variable(2).value(i);
            figure(1);
            set(gcf, 'Position', [125   319   549   470], 'Color', [1 1 1],...
                'PaperPositionMode', 'auto');
            %draw velocity
            results.flow.current = 2;
            drawNormalizedVelocity (gca, results, max_min_velocity, scale_velocity);
            title('Velocity field [mm.s^{-1}]', 'Fontsize', 12);
            xlabel(sprintf('time %0.1f day', t/24), 'Fontsize', 12);
            colormap gray;
            %draw shear
            results.flow.current = 1;
            figure(2);
            set(gcf, 'Position', [125   319   549   470], 'Color', [1 1 1],...
                'PaperPositionMode', 'auto');
            drawNormalizedVelocity (gca, results, max_min_shear, scale_shear);
            title('Shear stress [N.m^{-2}]', 'Fontsize', 12);
            xlabel(sprintf('time %0.1f day', t/24), 'Fontsize', 12);
            colormap gray;

            %save
            fn1 = sprintf('%s/frames/flow%06d.png', baseDir, i);
            saveFigure(1, fn1, [baseDir '/frames/']);
            fn2 = sprintf('%s/frames/shear%06d.png', baseDir, i);
            saveFigure(2, fn2, [baseDir '/frames/']);
            % join the two images
            fn3 = sprintf('%s/frames/frm%06d.png', baseDir, i);
            joinImages(fn1, fn2, fn3);
            delete(fn1);
            delete(fn2);
        catch
            disp(lasterr);
        end;
    end;
end;