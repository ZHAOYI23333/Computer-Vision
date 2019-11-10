% Author: Yi Zhao
% CSE 5524, HW1
% 09/01/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1

a = 0.4;
wx = [.25-.5*a .25 a .25 .25-.5*a];
wy = reshape(wx, [5,1]);
% double, grayscale image. Make sure it's size is (odd, odd)
gauss = double(rgb2gray(imread('./data/1.jpg'))); 
for i = 1:4
    % save gaussian pyramid from layer1 to layer2
    imwrite(gauss/255, sprintf('./output/pyramid/gaussian_%d.png',i));
    nextGaussian = blurNsample(gauss, wx, wy);
    est = interpolation(nextGaussian);
    imwrite((gauss - est)/255, sprintf('./output/pyramid/laplacian_%d.png',i));
    % Set gauss to nextGaussian to continue building pyramid
    gauss = nextGaussian;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2

T_list = [10 20 35 50 75 100 125 150 200 255];
backgroundIm = double(imread('./data/bg000.bmp'));
inputIm = double(imread('./data/walk.bmp'));
binary_output = abs(inputIm - backgroundIm);
for T = T_list
    binary_copy = binary_output;
    binary_copy(binary_copy >= T) = 1;
    binary_copy(binary_copy ~= 1) = 0;
%     imshow(binary_copy)
%     pause;
    imwrite(binary_copy, sprintf('./output/background_sub1/T%d.png',T))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3

% Read images. Store 30 images as a cube 
for i = 1:30
    filename = sprintf('./data/bg%03d.bmp', i-1);
    Im(:,:,i) = double(imread(filename));
end
% Calculate stdev & mean matrix
stdev = std(Im, 0, 3);
avg = mean(Im, 3);
% Do background subtraction 2
bn_output = (inputIm - avg).^2 ./ stdev.^2;
T_list2 = 1:20;
for T = T_list2
    bn_copy = bn_output;
    bn_copy(bn_copy >= T^2) = 1;
    bn_copy(bn_copy ~= 1) = 0;
%     imshow(bn_copy)
%     pause;
    imwrite(bn_copy, sprintf('./output/background_sub2/T%d.png',T))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4

% Based on observation, 'T250.0.png' is used for best output
bsIm = double(imread('./output/background_sub2/T16.png'));
d_bsIm = bwmorph(bsIm, 'dilate');
imwrite(d_bsIm, './output/dilated_ouput.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 5

[L, num] = bwlabel(d_bsIm, 8);
% select second most label (the most is 0)
second_most = mode(L(L ~= 0), 'all');
L(L ~= second_most) = 0;
imwrite(L, './output/final_ouput.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 Helper function

% First blur the image using gaussian separable filter wx and wy, 
% and then sample the image to half size.
function bluredSampled = blurNsample(image, wx, wy)
    blured = imfilter(imfilter(image, wy, 'replicate'), wx, 'replicate');
    bluredSampled = blured(1:2:end, 1:2:end);
end

% Implement the interpolation function
function newImage = interpolation(image)
    newImage = zeros(size(image)*2-1);
    newImage(1:2:end, 1:2:end) = image;
    row_ave = conv2(image, [1 1], 'valid')/2;
    col_ave = conv2(image, [1;1], 'valid')/2;
    newImage(1:2:end,2:2:end) = row_ave;
    newImage(2:2:end,1:2:end) = col_ave;
    mid_ave = conv2(row_ave,[1;1], 'valid')/2;
    newImage(2:2:end, 2:2:end) = mid_ave;
end

