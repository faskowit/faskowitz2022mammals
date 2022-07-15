function [ netsimile ] = netsimile_dist_bu(amat,bmat)
% weighted implementation of netsimile algo for comparing networks
% https://arxiv.org/abs/1209.2684

% enfore binary
amat = amat > 0 ;
bmat = bmat > 0 ; 

% get the netsimile of each network
ns_a = feat_extract(amat) ;
ns_b = feat_extract(bmat) ;

% canbera distance
canbdist = @(a_,b_) sum(abs(a_-b_)./(abs(a_)+abs(b_)),'omitnan') ;
netsimile = canbdist(ns_a,ns_b) ;

end 

function [ nsfeat ] = feat_extract(imat_) 

n = size(imat_,1) ;
shpath = distance_bin(imat_>0) ;
triunroll = @(x_) x_(logical(triu(ones(size(x_,1)),1))) ;

% seven features
featnames = { 'd' 'cc' 'd2' 'c1' 'ego_e' 'ego_o' 'ego_n' } ;
feat_str = struct() ;

% weighted degree
feat_str.d = degrees_und(imat_)' ;
% clusterin coeff
feat_str.cc = clustering_coef_bu(imat_) ;
% two-hop neighbors degree
feat_str.d2 = nan(n,1) ;
for ii = 1:n
    nei2 = shpath(ii,:)==2 ;
    feat_str.d2(ii) = mean(feat_str.d(nei2)) ;  
end
% average clust of neighborhood
feat_str.c1 = nan(n,1) ;
for ii = 1:n
    nei = shpath(ii,:)==1 ;
    feat_str.c1(ii) = mean(feat_str.cc(nei)) ;  
end

% egonet stuff
feat_str.ego_e = nan(n,1) ;
feat_str.ego_o = nan(n,1) ;
feat_str.ego_n = nan(n,1) ;

% loop through egonet for each node
for ii = 1:n
    nei = shpath(ii,:)==1 | shpath(ii,:)==0 ; % path == 0 to include self
    egon = imat_(nei,nei) ;
    
    % sum edges in egonet
    feat_str.ego_e(ii) = sum(triunroll(egon)) ;
    
    % sum edges outgoing egonet
    feat_str.ego_o(ii) = sum(imat_(nei,~nei),'all') ;
    
    % get the neighbors of the egonet
    nei2 = any(shpath(nei,:) == 1) ;
    % get rid of the edges within within-egonet connections
    nei2 = nei2 - nei ;
    
    % weighted degree of neighbors of egonet
    feat_str.ego_n(ii) = mean(feat_str.d(logical(nei2))) ;
end

% get features 
nsfeat = nan(35,1) ;
for ii = 0:(length(featnames)-1)
    jj = (ii*5)+1 ;
    kk = jj+4 ; 
    nsfeat(jj:kk) = feat_agg(feat_str.(featnames{ii+1})) ;
end

end

function ff = feat_agg(iarr)

ff = nan(5,1) ;
ff(1) = median(iarr) ;
ff(2) = mean(iarr) ;
ff(3) = std(iarr) ;
ff(4) = skewness(iarr) ;
ff(5) = kurtosis(iarr) ;

end
