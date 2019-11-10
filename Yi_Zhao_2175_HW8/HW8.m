% Author: Yi Zhao
% CSE 5524, HW8
% 10/26/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1

image = double(imread("./data/checker.png"));
sigma_window = 1; sigma_grad = 0.7; alpha = 0.05; T = 1000000;

% gaussian window with sigma = 1, s * sigma mask
G = fspecial('gaussian', 2*ceil(3*sigma_window)+1, sigma_window);
[Gx, Gy] = gaussDeriv2D(sigma_grad);

Ix = imfilter(image, Gx, 'replicate');
Iy = imfilter(image, Gy, 'replicate');
Ixy = Ix .* Iy;

gIx2 = imfilter(Ix.^2, G, 'replicate');
gIy2 = imfilter(Iy.^2, G, 'replicate');
gIxy = imfilter(Ixy, G, 'replicate');

R = gIx2.* gIy2 - gIxy.^2 - alpha*(gIx2+gIy2).^2;

disp(R(17:23, 17:23))

R(R<T) = 0;
imagesc(R)
colormap('gray')
title('Thresholded R','FontSize', 14)
pause;
saveas(gcf,'./output/ThresholdedR.jpg')

duplicate = R;
for r = 1:size(R, 1)
    for c = 1:size(R, 2)
        w = getWindow(3, duplicate, r, c);
        if duplicate(r, c) ~= max(w, [], 'all')
            R(r, c) = 0;
        end
    end
end

imshow(image/255, 'InitialMagnification','fit')
hold on
[y, x] = find(R);
plot(x, y, 'r.', 'MarkerSize', 20)
title('Harris Detector after Non-Maximum Suppression', 'FontSize',14)
saveas(gcf,'./output/Detector.jpg')
hold off
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2

img = double(imread("./data/tower.png"));
T_list = [10 20 30 50]; n = 9;

for T = T_list
    fast = getFAST(img, T, n);
    imshow(img/255, 'InitialMagnification','fit')
    hold on
    [y, x] = find(fast);
    plot(x, y, 'r.')
    title(sprintf('FAST on Tower with T = %i', T), 'FontSize',14)
    saveas(gcf, sprintf('./output/FAST_T%i.jpg', T))
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Function

% Gaussian derivative function from HW2
function [Gx, Gy] = gaussDeriv2D(sigma)
    length = 2 * ceil(sigma * 3) + 1;
    for r = 1:length
        for c = 1:length
            y = -r + ceil(3*sigma) + 1;
            x = c - ceil(3*sigma) - 1;
            Gx(r,c) = -x * exp(-1 * (x^2 + y^2)/(2 * sigma.^2)) / (2 * pi * sigma^4);
            Gy(r,c) = -y * exp(-1 * (x^2 + y^2)/(2 * sigma.^2)) / (2 * pi * sigma^4);
        end
    end
end


% Get the square surrounded img(r, c) with length = len
function window = getWindow(len, img, r, c)
    left = max(1, c - floor(len/2));
    right = min(size(img, 2), c + floor(len/2));
    top = max(1, r - floor(len/2));
    bottom = min(size(img, 1), r + floor(len/2));
    window = img(top:bottom, left:right);
end


% Hardcode the border with r = 3. 
function border = getBoarder(img, r, c)
    top = img(r-3, c-1:c+1);
    bottom = img(r+3, c-1:c+1);
    left = img(r-1:r+1, c-3);
    right = img(r-1:r+1, c+3);
    lt = img(r-2, c-2);
    rt = img(r-2, c+2);
    lb = img(r+2, c-2);
    rb = img(r+2, c+2);
    border = [top, rt, right', rb, flip(bottom), lb, flip(left'), lt];
end


% Check if the array satisfy n >= n*.
function res = verifyList(category, n)
    list = [category, category];
    l = 1; temp = 1;
    for i = 2:size(list, 2)
        if list(i-1) == list(i) && list(i) ~= 0
            temp = temp + 1;
        else
            if temp > l
                l = temp;
            end
            temp = 1;
        end
    end
    res = l >= n;
end

% Get FAST feature points matrix
function fast = getFAST(img, T, n)
    fast = zeros(size(img));
    for r = 4:size(img,1)-3
        for c = 4:size(img,2)-3
            b = getBoarder(img, r, c);
            category = zeros(size(b));
            category(b > img(r,c)+T) = 1;
            category(b < img(r,c)-T) = -1;
            fast(r,c) = verifyList(category, n);
        end
    end
end
 
