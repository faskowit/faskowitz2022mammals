function c = fcn_mst_plus(a,dens,dist)
% compute threhsold with minimum spanning tree and edges added to reach
% certain percent density
% inputs: matrix, density [0-1], distance transformation (inv, minus)
% outputs: binary mask to do thresholding
% written by Rick Betzel, Josh Faskowitz additions, Indiana University,
% 2021

if nargin < 3
    dist = 'inv' ;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch dist
    case 'inv'
        dist = 1./a ; % use for structual matrix
    case 'minus'
        dist = 1 - a ; % use for correlation matrix
    otherwise
        error([ 'dist option: ' dist ' not valid' ]) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amst = graphminspantree(sparse(dist),'method','kruskal');
amst = amst + amst';

n = length(a);
m = n*(n - 1)/2;
tgt = round(m*dens);
mmst = nnz(triu(amst));
add = tgt - mmst;
[u,v,w] = find(a.*triu(a & ~amst,1));
[~,idxsort] = sort(w,'descend');
idx = (v - 1)*n + u;
b = amst;
try
   b(idx(idxsort(1:add))) = 1;
catch
   error([ 'not enough edges to reach desired density of: ' num2str(dens) ])
end
c = double(b | b');

