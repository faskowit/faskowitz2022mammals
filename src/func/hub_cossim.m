function [ cosdiff, origcos ] = hub_cossim(mat,distmat,hubmethod,Nreps)
% the probability that hubs remain after 100% rewire

if nargin < 4
   Nreps = 100 ; 
end

% n = size(mat,1) ;
% enfore und
mat = (mat + mat') ./ 2 ;

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% original hubs
[~,orighub] = get_hub_score_wei_und(mat,hubmethod) ;

o_cos = squareform(1-pdist(mat,'cosine')) ;
o_hubcos = mean(triunroll(o_cos(orighub,orighub)),'all') ;
o_feedcos = mean(o_cos(~orighub,orighub),'all') ;
o_othcos = mean(triunroll(o_cos(~orighub,~orighub)),'all') ;

rw_hubcos = nan(Nreps,1) ;
rw_feedcos = nan(Nreps,1) ;
rw_othcos = nan(Nreps,1) ;

for ndx = 1:Nreps
    disp(ndx)
    % full (100%) rewire using Gs
    [~,rewr_mat] = geombinsurr_partial(mat,distmat,1,20,'quantiles') ;
    % find the hubs
    [~,hh] = get_hub_score_wei_und(rewr_mat,hubmethod) ;

    rep_cos = squareform(1-pdist(rewr_mat,'cosine')) ;
    
    rw_hubcos(ndx) = mean(triunroll(rep_cos(hh,hh)),'all') ;
    rw_feedcos(ndx) = mean(rep_cos(~hh,hh),'all') ;
    rw_othcos(ndx) = mean(triunroll(rep_cos(~hh,~hh)),'all') ;
    
end

cosdiff = nan(3,1) ;
cosdiff(1) = o_hubcos - median(rw_hubcos) ;
cosdiff(2) = o_feedcos - median(rw_feedcos) ;
cosdiff(3) = o_othcos - median(rw_othcos) ;

origcos = nan(3,1) ;
cosdiff(1) = o_hubcos ;
cosdiff(2) = o_feedcos  ;
cosdiff(3) = o_othcos ;
