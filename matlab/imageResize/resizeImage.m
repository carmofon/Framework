function [] = resizeImage(imFileName, f)

% [] = resizeImage(imFileName, f)
% resizes and image. Srhinks image if f < 1 and enlarges if f >1
% creates a new file of same format
% and appends to file name "_big" or "_small" 


if (f > 1)
    % if n > 1 appends '_big' to filename
    add = '_big';
elseif (f < 1)
    % if n < 1 appends '_small' to filename
    add = '_small';
else
    % if n ==1 give error
    error('Trying to rescale image with f = 1');
end


% load image
im = imread(imFileName); 

% decompose filename to get the image format
% assumes that last 3 characters in filename are the file's extension
type = imFileName(end-2:end);
name = imFileName(1:(end-4));

% get sizes
[m, n, k] = size(im); 
m1 = round(m*f);
n1 = round(n*f);


% iterates trhough each of the RGB color channels and resizes each
imSmall = [];
for i = 1:k
    [X,Y] = meshgrid( (0:n-1)/(n-1), (0:m-1)/(m-1) );
    % new sampling location
    [XI,YI] = meshgrid( (0:n1-1)/(n1-1) , (0:m1-1)/(m1-1) );
    imSmall(:, :, i) = interp2( X, Y, im(:, :, i), XI,YI ,'cubic'); % the new image
end;

% reconvert image from double to 8-bit per channel
imSmall = uint8(imSmall);
 
% write file
imwrite(imSmall, [name add '.' type], type);