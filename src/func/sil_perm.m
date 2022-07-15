function [sil,pval] = sil_perm(points,ca,numP)

% the empirical eval
tmp = evalclusters(points,ca,meth) ;
sil = tmp.CriterionValues ;

% scramble 
% ca_shuff = zeros(length(ca),numP) ;
res = zeros(numP,1) ;
shuff = @(v_) v_(randperm(length(v_))) ;
for idx = 1:numP
    if mod(idx,100)==0 ; disp(idx) ; end
    tmp = evalclusters(points,shuff(ca),meth)  ;
    res(idx) = tmp.CriterionValues ;
end

pval = (sum(res>=sil)+1) / ( numP + 1 ) ;
