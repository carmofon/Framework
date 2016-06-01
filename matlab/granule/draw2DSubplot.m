function [] = draw2DSubplot(data2d, figureNo, subplotNo)

figure(figureNo);
set(figureNo, 'Position', [1          31        1280         918]);
%subplot(4, 6, subplotNo);
try,
    imagesc(data2d);
    %contourf(data2d);
    colormap gray;
    colorbar;
    axis equal tight;
    set(gca, 'XTick', [], 'YTick', []);
catch
end;