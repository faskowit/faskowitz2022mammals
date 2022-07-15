% function fcn_gnets(
clear all
close all
clc

k = 5;
params = rand(k*(k - 1),1)*10;
N = 128;
sz = rand(k,1);
sz = sz/sum(sz);
sz = round(sz*N);
while sum(sz) ~= N
    dff = sum(sz) - N;
    if dff > 0
        for i = 1:dff
            r = randi(k);
            if sz(r) > 0
                sz(randi(k)) = sz(randi(k)) - 1;
            end
        end
    else
        for i = 1:-dff
            r = randi(k);
            sz(randi(k)) = sz(randi(k)) + 1;
        end
    end
end

ci = zeros(N,1);
count = 0;
for j = 1:k
    idx = (count + 1):(count + sz(j));
    ci(idx) = j;
    count = count + sz(j);
end

omega = triu(rand(k));
omega = omega + triu(omega,1)';
m = 1000;
P = omega(ci,ci);
mask = triu(ones(N),1) > 0;
Pmask = P(mask);
Pmask = m*Pmask./sum(Pmask);
Q = zeros(N);
Q(mask) = Pmask;
Q = Q + Q';
omega = zeros(k);
for i = 1:k
    for j = 1:k
        vals = Q(ci == i,ci == j);
        omega(i,j) = unique(nonzeros(vals));
    end
end

A = zeros(N);
count = 0;
for i = 1:k
    for j = i:k
        count = count + 1;
        if i == j
            mask = rand(sz(i));
            a = triu(mask <= omega(i,j),1);
        else
            mask = rand(sz(i),sz(j));
            a = mask <= omega(i,j);
        end
        idx = find(a);
        m = length(idx);
        w = exprnd(params(count,1),[m,1]);
        g = zeros(size(a));
        g(idx) = w;
        A(ci == i,ci == j) = g;
    end
end
A = A + A';