function [coef] = myassort_und(mat,corrtype)

if nargin<2
   corrtype = 'pearson' ;
end

tmpdist = sum(mat,2,'omitnan');
[i,j] = find(triu(mat,1)>0);
ii = tmpdist(i);
jj = tmpdist(j);

coef = corr(ii,jj,'type',corrtype) ;
