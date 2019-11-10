% Author: Yi Zhao, zhao.2175@osu.edu
% Professor: Jim Davis
% Class: CSE 5524, MW 12:45PM
% Date: 09/06/2019
%%%%%%%%%%%%%%%%%%%% HW 2 %%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1

sigma=20;  % use different values
G = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
faceIm=double(imread('./data/affleck_gray.png'));
gIm = imfilter(faceIm, G, 'replicate');
imshow(gIm/255);  % double images need range of 0-1
imwrite(uint8(gIm), './output/gIm.bmp');
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 (Function definition locates at the end)

[a,b] = gaussDeriv2D(5);  % change sigma when needed
img_gx = imagesc(a);
colormap('gray');
pause;

img_gy = imagesc(b);
colormap('gray');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3

[Gx,Gy] = gaussDeriv2D(5);
Im = double(imread('./data/img.jpg'));
gxIm = imfilter(Im, Gx, 'replicate');
gyIm = imfilter(Im, Gy, 'replicate');
magIm = sqrt(gxIm.^2 + gyIm.^2);
imagesc(magIm);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4
T = 3;
tIm = magIm > T;
imagesc(tIm);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 5

T = 100;
Fx = -fspecial('sobel')';
fxIm = imfilter(Im,Fx);
Fy = -fspecial('sobel');
fyIm = imfilter(Im,Fy);
magIm = sqrt(fxIm.^2 + fyIm.^2);
tIm = magIm > T;
imagesc(tIm);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 6
Im = rgb2gray(Im);
edge(Im, 'canny')
pause;
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Function definition 

function [Gx, Gy] = gaussDeriv2D(sigma)
    length = 2 * ceil(sigma * 2) + 1;
    for r = 1:length
        for c = 1:length
            y = -r + ceil(2*sigma) + 1;
            x = c - ceil(2*sigma) - 1;
            Gx(r,c) = -x * exp(-1 * (x^2 + y^2)/(2 * sigma.^2)) / (2 * pi * sigma^4);
            Gy(r,c) = -y * exp(-1 * (x^2 + y^2)/(2 * sigma.^2)) / (2 * pi * sigma^4);
        end
    end
end


