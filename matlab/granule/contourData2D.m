function [] = contourData2D(data2d, label)

contourf(data2d(end:-1:1, :));
xlabel(label);
colormap gray, colorbar('SouthOutside');
axis image;
set(gca, 'XTick', [], 'YTick', []); 