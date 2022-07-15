function comb = f_NcombK(N,K,r)

comb = zeros(r,K);

for i=1:r
    rr = randperm(N);
    comb(i,:) = rr(1:K);
end