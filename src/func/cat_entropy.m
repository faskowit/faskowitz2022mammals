function [e,h] = cat_entropy(ci)
if ~all(floor(ci) == ci)
    error('need category data')
end
ele = unique(ci)' ;
countcount = @(AA,BB) arrayfun(@(x)sum(AA == x), BB) ;
h = countcount(ci,ele) ;
p = h ./ sum(h,2);
e = -sum(p.*log2(p),2,'omitnan');
