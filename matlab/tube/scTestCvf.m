% scTestCvf - routines for debbuging the cvf

iterationNumber = 40;
baseDir = 'C:\Documents and Settings\xavier\Desktop\results\tube';

% read the map
fn = sprintf('%s/cvf/map%08d.txt', baseDir, iterationNumber);
map =  load(fn, 'ascii');

% read the cvf
fn = sprintf('%s/cvf/cvf%08d.txt', baseDir, iterationNumber);
cvf =  load(fn, 'ascii');

% read the init
fn = sprintf('%s/cvf/init%08d.txt', baseDir, iterationNumber);
init =  load(fn, 'ascii');

figure(1);
set(gcf, 'Position', [148 50 1060 840], 'Color', [1 1 1],...
    'PaperPositionMode', 'auto');
%draw MAP
imagesc(map + 100*cvf);
title('map and CVF');
colorbar;
axis image;
