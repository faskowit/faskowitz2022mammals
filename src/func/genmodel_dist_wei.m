function [maxD,meanD,K] = genmodel_dist_wei(aij,bij,d_aij,d_bij)
% the distance function from generative modeling paper

% note the binary functions
evalfuncs = cell(4,1) ;
evalfuncs{1} = @(a_) sum(a_,2) ;
evalfuncs{2} = @(a_) clustering_coef_wu(a_);
evalfuncs{3} = @(a_) my_btwn_cent(a_)';
evalfuncs{4} = @(a_,d_) d_(triu(a_,1) > 0);

whicheval = 1:4 ; 
nEval = length(whicheval) ; 
K = zeros(nEval,1);

for ee = 1:nEval
    
    if ee < 4
        x = evalfuncs{ee}(aij) ;
        y = evalfuncs{ee}(bij) ;
    else % do the distance one
        x = evalfuncs{ee}(aij,d_aij) ;
        y = evalfuncs{ee}(bij,d_bij) ;
    end
    
    K(ee) = fcn_ks(x(:),y(:));
end
    
maxD = max(K);
meanD = mean(K);

end

function kstat = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);

end

% added better betweenness centrality
function btwncent = my_btwn_cent(W,transdist)
% using floyd. do the inverse weight inside func

if nargin < 2
    transdist = 'inv' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch transdist
    case 'inv'
        L = 1 ./ W;
    case 'log'
        L = -log(W./ (max(W)+eps)) ; % +eps makes no 0 distances
    otherwise
        error('not a valid transdist') 
end

n = size(W,1) ;
[~,hops_,paths_] = distance_wei_floyd(L) ;
btwncent = zeros(n,1) ; 

for idx = 1:n
    for jdx = 1:n
        if idx >= jdx ; continue ; end
        
        pa = retrieve_shortest_path(idx,jdx,hops_,paths_) ;
        btwncent(pa(2:(end-1))) = btwncent(pa(2:(end-1)))+2 ;
    end
end

end % end of my betweenness
