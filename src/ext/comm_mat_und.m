function [B,Bnorm,commsizes,blockpairs] = comm_mat_und(W,ci,omitnan)
%COMM_MAT 		Returns the block matrix of the community structure of a graph W with given membership vector ci
%
%   Inputs      W,  undirected weighted/binary network. 
%               ci, membership vector
%
%   Outputs:    B,  block matrix of community structure. Element B(i,j) contains the sum of links 
%				from community to community j, half of them when i==j (no self-loops are considered).
%               Bnorm, bloc matrix normalized by size of blocks (i.e. the
%               numbers in blockpairs)
%               commsizes, size of communities
%               blockpairs, number of possible edges between blocks
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%  https://github.com/CarloNicolini/communityalg/blob/master/comm_mat.m
% joshedit

if any(diag(W))
	warning('Input matrix has self-loops');
end

if nargin < 3
   omitnan = 1 ; 
end

ci=ci(:)'; % force it to be a row vector
if omitnan % get rid of NaN?
   W(isnan(W)) = 0 ;
end

% The first matrix C is an ncommsX lenght(ci) matrix of ones and zeros. 
% Each row holds a 1's where ci is equal to one of the unique value.
C=bsxfun(@eq,ci,unique(ci)');
% Obtain the block matrix
B = C*W*C';
% Divide the diagonal because self-loops are counted twice
B(logical(eye(size(B)))) = B(logical(eye(size(B))))./2;

% Compute also normalized version of block matrix
commsizes = sum(C,2);
% number of pairs inside every community (no self loops) on diagonal
% ni*(ni-1)/2
commpairs = diag(commsizes)*diag(commsizes-1)/2;
% number of pairs between pairs of communities ni*nj
commpairs2=diag(commsizes)*ones(length(commsizes))*diag(commsizes);
blockpairs = commpairs2.*(1-eye(length(commsizes))) + commpairs;
% on the diagonal the intracommunity densities, outside the diagonal the
% intercommunity densities
Bnorm = B./blockpairs;

% n= length(commsizes) ;
% sanity = zeros(n) ;
% comms = unique(ci) ;
% for idx = 1:n
%     for jdx = 1:n
%         if idx > jdx ; continue 
%         elseif idx == jdx 
%             sanity(idx,jdx) = (sum(W(ci==comms(idx),ci==comms(jdx)),'all','omitnan'))/2 ;
%         else
%             sanity(idx,jdx) = sum(W(ci==comms(idx),ci==comms(jdx)),'all','omitnan') ;
%         end
%     end
% end
% sanity = sanity + triu(sanity,1)' ;

