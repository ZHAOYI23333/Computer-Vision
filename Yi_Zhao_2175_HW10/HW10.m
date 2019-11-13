% Author: Yi Zhao
% CSE 5524, HW10
% 11/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1

left = double(imread('./data/left.png'));
right = double(imread('./data/right.png'));
map = zeros(size(left));
[r,c] = size(left);
l = 11; % window size
h = floor(l/2);

for y = 1:r-l+1
    for x = 1:c-l+1
        leftwindow = left(y:y+11-1, x:x+11-1);
        ncclist = zeros(1,x-max(1, x-50)+1);
        for z = 1:x-max(1, x-50)+1
            rightwindow = right(y:y+11-1, x-z+1:x-z+11);
            ncclist(z) = calculateNCC(leftwindow, rightwindow);
        end
        [m, m_ind] = max(ncclist);
        map(y+h, x+h) = m_ind;
    end
end
figure;
imagesc(map, [0 50]);
axis equal;
colormap("gray")
title("The Calculated Disparse Map")
saveas(gcf,"./output/depth.png");
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2 KNN
load('data/train.txt')
load('data/test.txt')
x_train = train(:,1:2);
y_train = train(:,3);
x_test = test(:, 1:2);
y_test = test(:, 3);

K = [1 5 11 15];

for k = K
   predict = KNN(x_train, y_train, x_test, y_test, k);
   plot(x_test(predict==y_test & predict == 1,1), x_test(predict==y_test & predict == 1,2), 'b.')
   axis equal;
   hold on
   plot(x_test(predict==y_test & predict == 2,1), x_test(predict==y_test & predict == 2,2), 'r.')
   plot(x_test(predict~=y_test,1 ), x_test(predict ~= y_test,2), 'ko')
   hold off
   title(sprintf('KNN Result on K = %d',k),'FontSize',14)
   saveas(gcf,sprintf("./output/k%d.png",k));
   pause;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Function

% return the value of NCC of two eqo_sized image.
function ncc = calculateNCC(origin, template)
    p_mean = mean(origin, [1 2]);
    t_mean = mean(template, [1 2]);
    p_sigma = std(origin, 0, [1 2]);
    t_sigma = std(template, 0, [1 2]);
    arr = ((origin - p_mean).*(template - t_mean))./(p_sigma.*t_sigma);
    % prevent inf value from dividing 0
    arr(p_sigma * t_sigma == 0) = min(arr, [], 'all'); 
    ncc = sum(sum(arr, [1 2])/(size(template,1)*size(template,2)-1), 'all');
end

% return a column of predicted labels
function predict = KNN(x_train, y_train, x_test, y_test, k) 
    index = knnsearch(x_train, x_test, 'K', k);
    predict = y_train(index);
    predict = mode(predict, 2);
    accuracy = sum(predict == y_test)/size(y_test, 1);
    fprintf('The accuracy at K = %u is: %.3f%%\n', k, accuracy*100)
end

