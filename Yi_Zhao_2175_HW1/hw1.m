% Author: Yi Zhao
% CSE 5524, HW1
% 09/01/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1

grayIm = imread('./data/buckeyes_gray.bmp');
imagesc(grayIm);
axis('image');
colormap('gray');
imwrite(grayIm, './output/buckeyes_gray.jpg'); 
pause;

rgbIm = imread('./data/buckeyes_rgb.bmp');
imagesc(rgbIm);
axis('image');
imwrite(rgbIm, './output/buckeyes_rgb.jpg');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2

grayIm = rgb2gray(rgbIm);
imagesc(grayIm);
axis('image');
imwrite(grayIm, './output/buckeyes_gray_converted.jpg'); 
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3
zBlock = zeros(10,10);
oBlock = ones(10,10)*255;
pattern = [zBlock oBlock; oBlock zBlock]; 
checkerIm = repmat(pattern, 5, 5); 
imwrite(uint8(checkerIm), './output/checkerIm.bmp'); 
Im = imread('./output/checkerIm.bmp');
imagesc(Im)
colormap('gray')
axis('image');
