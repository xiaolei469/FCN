function p = pX(numx,numy,N,sX)
% Kirsch A, Mitzenmacher M, Pietracaprina A, et al. An efficient rigorous approach for identifying statistically significant frequent itemsets[J].
% Journal of the ACM (JACM), 2012, 59(3): 1-22.
p = 0;
for n=sX:N
    p = p + pXn(numx,numy,N,n);
end
end

function pn = pXn(numx,numy,N,n)
log_pn = nchoosekln(N,n) + n*(log(numx*numy)-2*log(N)) + (N-n)*(log(N^2-numx*numy)-2*log(N));
pn = exp(log_pn);
end

function nk = nchoosekln(n,k)
nk = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
end



