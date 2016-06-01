% scDrawFig1_phototrophic

dirName = 'E:\jxavier\sloughingPaper\data';
dirOut = 'E:\jxavier\sloughingPaper\figs_temp/';
fileName = '/r01l1i11days.txt';

growthCurveData = readGrowthCurve([dirName fileName]);


figure(1);
set(1, 'Position', [40 420 1202 515], 'PaperPositionMode', 'auto');
for i = 1:9,
    subplot(3, 6, floor((i-1)/3)*6 + mod(i-1,3) + 1);
    plot(growthCurveData.time, growthCurveData.sensor(i).value, 'k-');
    xlabel('time [day]');
    ylabel('X');
end;

% PHASE-SPACE
% group the data into [X_{N+1}; X_{N}] pairs
subplot(1, 2, 2);
pairs = [];
for i = 1:length(growthCurveData.sensor),
    x1 = growthCurveData.sensor(i).value(1:end-1);
    x2 = growthCurveData.sensor(i).value(2:end);
    if isempty(pairs),
        pairs = [x1'; x2'];
    else
        pairs = [pairs(1,:), x1'; pairs(2,:), x2'];
    end;
    h = plot(x1, x2, 'k-');
    hold on;
end;
hold off;
xlabel('X_n');
ylabel('X_{n+1}');

print('-depsc', [dirOut 'experiment.eps']); 
