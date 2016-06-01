function [] = joinImages(fn1, fn2, fn3)

im1 = imread(fn1);
im2 = imread(fn2);

[m, n, q] = size(im1);

%im3 = zeros(m, 2*n, q);
im3(:,1:n,:) = im1;
im3(:,n+1:2*n,:) = im2;

imwrite(im3, fn3, 'png');