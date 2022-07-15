function [weiBM,avgWeiBM,binBM,avgBinBM,sizesMat] = get_block_mat(cij,ca,noNan)
% given an adjacency matrix + community affiliations, return a block matrix
%
% inputs: 
%           cij:        adjacency matrix
%           ca:         community affiliation vector
%           noNan:      do not count nan edges when getting size of blocks
%
%           NOTE: this function will treat NaNs as zeros
%
% returns: 
%           weiBM:      sum of weights block matrix
%           avgWeiBM:   average weighted block matrix
%           binBM:      sum of binary weights block matrix
%           avgBinBM:   average binary block matrix
%           stdWeiBM:   std weights block matrix 
%           sizesMat:   number of edges per block (incl. NaN edges)
%
% Josh Faskowtiz IU

% make ca column vec
if ~iscolumn(ca)
   ca = ca'; 
end

if nargin < 3
   noNan = 0 ; 
end

% number coms
nBlocks = length(unique(ca));

% number nodes per block
% blockSizes = histcounts(sort(ca));
countcount = @(v_) arrayfun(@(x_)sum(v_(:) == x_), unique(v_)) ;
blockSizes = countcount(ca) ;

sizesMat = bsxfun(@times,...
    (bsxfun(@times,ones(nBlocks),blockSizes)),...
    (blockSizes)');

% % assume we ignore self connections
sizesMat(~~eye(nBlocks)) = sizesMat(~~eye(nBlocks)) - blockSizes(:) ;

W = cij ;

% get rid of nans
W(isnan(W)) = 0 ;

C = dummyvar(ca) ;
C = C' ;
% C2 = bsxfun(@eq,ca,unique(ca)');
B = C*W*C' ;
Bcounts = C*(W>0)*C' ;

if noNan == 1
   nanCounts = C*isnan(cij)*C'; 
   sizesMat = sizesMat - nanCounts ;
end

% weight block matrix
weiBM = B ;
% for the on the diagonal, we should not double count connections
% edit-> user should do this outside func
% weiBM(~~eye(nBlocks)) = weiBM(~~eye(nBlocks)) ./ 2 ;

% avg weight block matrix
avgWeiBM = B ./ sizesMat ;

% weight block matrix
binBM = Bcounts ;
% for the on the diagonal, we should not double count connections
% edit-> user should do this outside func
% binBM(~~eye(nBlocks)) = binBM(~~eye(nBlocks)) ./ 2 ;

% avg weight block matrix
avgBinBM = Bcounts ./ sizesMat ;