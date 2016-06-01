% scPlotVelocity

iteration = 19;
baseDir = 'C:\Documents and Settings\xavier\Desktop\results\tube\velocity\';
radiusTube = 1000;
gridSide = 65;
systemSize = radiusTube/(0.5 - 1 / gridSide) * 1.01;
voxelSide = systemSize/gridSide;



fileName=sprintf('%siteration%08d\\Velocity5.txt', baseDir, iteration);
v = load(fileName);

%compute the gradient
[stressX,stressY] = gradient(-v, voxelSide);

figure(1);

subplot(2, 2, 1);
contourf(v);
colorbar;
title('velocity [um/h]');
hold on
quiver(stressX, stressY, 2);
hold off 

subplot(2, 2, 2);
imagesc(sqrt(stressX.^2 + stressY.^2));
colorbar;
title('stressY [um/h]');

