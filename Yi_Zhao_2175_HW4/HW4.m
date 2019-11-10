% Author: Yi Zhao
% CSE 5524, HW4
% 09/21/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1

for i = 1:4
    % read image (in double): boxIm[1-4].bmp
    im = double(imread(sprintf('./data/boxIm%d.bmp',i)));
    % display 7 similituyde moment for each of image
    disp(similitudeMoments(im))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2

% Load the data
clear; close all;
load('./data/eigdata.txt');
X = eigdata;
subplot(2,1,1);
plot(X(:,1),X(:,2),'b.');
title('Data w/o Subtracted Mean','FontSize', 14)
axis('equal');
% mean-subtract data
m = mean(X);
Y = X - ones(size(X,1),1)*m;
subplot(2,1,2);
plot(Y(:,1),Y(:,2),'r.');
title('Data w/ Subtracted Mean','FontSize', 14)
axis('equal');
pause;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3

K = cov(Y);
[coeff,score,latent,~,explained] = pca(Y);
[U, V] = eig(K);

c = 9;
len1 = sqrt(c*V(1,1));
len2 = sqrt(c*V(2,2));

hold on
plot([U(1,1),0]*len1,[U(2,1),0]*len1, 'g','LineWidth',2)
plot([U(1,2),0]*len2,[U(2,2),0]*len2, 'r','LineWidth',2)
plot(Y(:,1),Y(:,2),'.')
title('Correlated Vector w/ 3 Std.','FontSize', 14)
ellipse(len1,len2,atan(U(1,1)/U(2,1)),0,0)
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 4

close all;

Y2 = Y * U';
plot(Y2(:,1),Y2(:,2),'.')
axis([-20 20 -20 20])
title('Unorrelated Y','FontSize', 14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1 Helper function

function Nvals = similitudeMoments(im)
    Nvals = [];
    
    % initialize matrix for row index, col index, x average and y average.
    xind = repmat(1:size(im,2),size(im,1),1); % col => x
    yind = repmat((1:size(im,1))', 1, size(im,2)); % row => y
   
    m00 = sum(im, 'all');
    m10 = sum(xind.*im, 'all');
    m01 = sum(yind.*im, 'all');

    xbar = ones(size(im)) * m10/m00;
    
    ybar = ones(size(im)) * m01/m00;
    % iteratively calculate 7 similitude moments
    for i = 0:3
        for j = max(0,(2-i)):(3-i)
            % 2 <= (i+j) <= 3
            nij = sum(((xind - xbar).^i).*((yind - ybar).^j).*im, 'all')/(m00.^((i+j)/2+1));
%             img = ((xind - xbar).^i).*((yind - ybar).^j).*im/(m00.^((i+j)/2+1));
%             imagesc(img)
%             colormap('gray')

%             title(sprintf('Similitude Moments: SM%d%d',i,j), 'FontSize', 14);
%             pause;
            Nvals = [Nvals, nij];
        end
    end
    
end
