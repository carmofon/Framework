function [] = saveFigure(fig, fn, tmpdir)

%save
figure(fig);
imFile = [tmpdir 'tmp.tif'];
print('-dtiff', '-r300', imFile);
% open image and resize it
im = imread(imFile);

MROWS = 293;
MCOLS = 343;

imSmall = imresize(im, [MROWS MCOLS], 'bilinear');
imwrite(imSmall, fn, 'png');
delete(imFile);
%
close(fig);
