% function pv = RunPoisson(X)
% % alpha = 0.2;
% [N,M] = size(X);
% Xb = Xbinaryeq2(X);
% sum_p = 0;
% num_alpha = 0;
% for m1=1:M-1
%     x1 = Xb{m1,1};
%     for m2=m1+1:M
%         x2 = Xb{m2,1};
%         x_1 = x1;
%         numx = length(find(x_1==1));
%         x_2 = x2;
%         numy = length(find(x_2==1));
%         % sX = (alpha*(numx+numy))/(1+alpha);
%         sX = 5;
%         sX = ceil(sX);
%         p = pX(numx,numy,N,sX);
%         inter = length(find(x_1.*x_2==1));
%         % ind =  double((inter)/(numx+numy-inter)>=alpha);
%         ind = double(inter >= sX);
%         sum_p = sum_p + sum(sum(p));
%         num_alpha = num_alpha + sum(sum(ind));
%     end
% end
% lambda = sum_p;
% pv = pvalue_poisson(num_alpha,lambda);
% end
% 
% function Xmb = Xbinaryeq2(X)
% [~,M] = size(X);
% Xmb = cell(M,1);
% for m=1:M
%     Xm = X(:,m);
%     Xb = double(Xm==2);
%     Xmb{m,1} = Xb;
% end
% end


% function pv = RunPoisson(X)
% [N,M] = size(X);
% 
% % Compute the degree of each node
% degree = sum(X == 2, 2);  % Sum along rows
% 
% % Compute the mean degree
% mean_degree = mean(degree);
% 
% % Set sX based on mean degree
% sX = 0.8 * mean_degree;
% sX = ceil(sX);
% 
% Xb = Xbinaryeq2(X);
% sum_p = 0;
% num_alpha = 0;
% 
% for m1=1:M-1
%     x1 = Xb{m1,1};
%     for m2=m1+1:M
%         x2 = Xb{m2,1};
%         x_1 = x1;
%         numx = length(find(x_1==1));
%         x_2 = x2;
%         numy = length(find(x_2==1));
% 
%         p = pX(numx,numy,N,sX);
%         inter = length(find(x_1.*x_2==1));
%         ind = double(inter >= sX);
% 
%         sum_p = sum_p + sum(sum(p));
%         num_alpha = num_alpha + sum(sum(ind));
%     end
% end
% 
% lambda = sum_p;
% pv = pvalue_poisson(num_alpha,lambda);
% end
% function Xmb = Xbinaryeq2(X)
% [~,M] = size(X);
% Xmb = cell(M,1);
% for m=1:M
%     Xm = X(:,m);
%     Xb = double(Xm==2);
%     Xmb{m,1} = Xb;
% end
% 
% end


% the code of k = 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pv = RunPoisson(X)
[N, M] = size(X);

% Compute common neighbors for all pairs
common_neighbors = zeros(M, M);
for m1=1:M-1
    for m2=m1+1:M
        common_neighbors(m1, m2) = sum((X(:, m1) == 2) & (X(:, m2) == 2));
    end
end

% Extract the upper triangular part (as the matrix is symmetric) and sort
sorted_common_neighbors = sort(nonzeros(triu(common_neighbors)), 'descend');

% Set sX based on top 0.1% of sorted_common_neighbors
index_1_percent = ceil(0.001 * length(sorted_common_neighbors));
sX = sorted_common_neighbors(index_1_percent);
disp([' sX ', num2str(sX)]);
Xb = Xbinaryeq2(X);
sum_p = 0;
num_alpha = 0;

for m1=1:M-1
    x1 = Xb{m1,1};
    for m2=m1+1:M
        x2 = Xb{m2,1};
        numx = sum(x1==1);
        numy = sum(x2==1);

        p = pX(numx,numy,N,sX);
        inter = sum(x1 .* x2 == 1);
        ind = double(inter >= sX);

        sum_p = sum_p + sum(p);
        num_alpha = num_alpha + sum(ind);
    end
end
lambda = sum_p;
disp(['lambda:', num2str(lambda)]);
disp(['num_alpha:', num2str(num_alpha)]);
pv = pvalue_poisson(num_alpha,lambda);
end

function Xmb = Xbinaryeq2(X)
[~, M] = size(X);
Xmb = cell(M,1);
for m=1:M
    Xmb{m,1} = double(X(:,m) == 2);
end
end


%{
% the code of k = 3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pv = RunPoisson(X)
    [N, M] = size(X);
    common_neighbors = zeros(M, M, M);
    for m1 = 1:M-2
        for m2 = m1+1:M-1
            for m3 = m2+1:M
                common_neighbors(m1, m2, m3) = sum((X(:, m1) == 2) & (X(:, m2) == 2) & (X(:, m3) == 2));
            end
        end
    end
    % sorted_common_neighbors = sort(nonzeros(triu(common_neighbors(:))), 'descend');
    % index_1_percent = ceil(0.1 * length(sorted_common_neighbors));
    % sX = sorted_common_neighbors(index_1_percent);
    sorted_common_neighbors = sort(nonzeros(common_neighbors(:)), 'descend');
    % disp(sorted_common_neighbors(1:10));
    if isempty(sorted_common_neighbors)
        disp('没有找到共有邻居，设置默认阈值 sX');
        sX = 0;  % 如果没有共有邻居或数据太少，可以设置一个默认阈值或处理为特殊情况
    else
        index_1_percent = ceil(0.00001 * length(sorted_common_neighbors));
        if index_1_percent == 0
            index_1_percent = 1;  % 确保索引至少为1
        end
        if index_1_percent > length(sorted_common_neighbors)
            index_1_percent = length(sorted_common_neighbors);  % 确保索引不超出数组长度
        end
        sX = sorted_common_neighbors(index_1_percent);
    end

    disp([' sX ', num2str(sX)]);
    Xb = Xbinaryeq2(X);
    sum_p = 0;
    num_alpha = 0;

    for m1 = 1:M-2
        x1 = Xb{m1, 1};
        for m2 = m1+1:M-1
            x2 = Xb{m2, 1};
            for m3 = m2+1:M
                x3 = Xb{m3, 1};
                numx = sum(x1 == 1);
                numy = sum(x2 == 1);
                numz = sum(x3 == 1);

                p = pX(numx, numy, numz, N, sX);
                inter = sum(x1 .* x2 .* x3 == 1);
                ind = double(inter >= sX);

                % sum_p = sum_p + p;
                % num_alpha = num_alpha + ind;
                % 与下面的写法实际上相同

                sum_p = sum_p + sum(p);
                num_alpha = num_alpha + sum(ind);
            end
        end
    end
    lambda = sum_p;
    disp(['lambda:', num2str(lambda)]);
    disp(['num_alpha:', num2str(num_alpha)]);
    pv = pvalue_poisson(num_alpha, lambda);
end

function Xmb = Xbinaryeq2(X)
[~, M] = size(X);
Xmb = cell(M,1);
for m=1:M
    Xmb{m,1} = double(X(:,m) == 2);
end
end
%}


% function pv = RunPoisson(X, sX)
% [N,M] = size(X);
% Xb = Xbinaryeq2(X);
% sum_p = 0;
% num_alpha = 0;
% for m1=1:M-1
%     x1 = Xb{m1,1};
%     for m2=m1+1:M
%         x2 = Xb{m2,1};
%         x_1 = x1;
%         numx = length(find(x_1==1));
%         x_2 = x2;
%         numy = length(find(x_2==1));
% 
%         p = pX(numx,numy,N,sX);
%         inter = length(find(x_1.*x_2==1));
%         % ind =  double((inter)/(numx+numy-inter)>=alpha);
%         ind = double(inter >= sX);
%         sum_p = sum_p + sum(sum(p));
%         num_alpha = num_alpha + sum(sum(ind));
%     end
% end
% lambda = sum_p;
% pv = pvalue_poisson(num_alpha,lambda);
% end
% 
% function Xmb = Xbinaryeq2(X)
% [~,M] = size(X);
% Xmb = cell(M,1);
% for m=1:M
%     Xm = X(:,m);
%     Xb = double(Xm==2);
%     Xmb{m,1} = Xb;
% end
% end
