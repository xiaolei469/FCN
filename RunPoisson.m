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
