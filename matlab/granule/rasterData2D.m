function [] = rasterData2D(data2d, label)

imagesc(data2d);
xlabel(label);
colormap gray, colorbar('SouthOutside');
axis image;
set(gca, 'XTick', [], 'YTick', []); 