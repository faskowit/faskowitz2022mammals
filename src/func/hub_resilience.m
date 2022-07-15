function [ hubresil, hubresilarea, resiljac, resilvar, orighub ] = ...
    hub_resilience(mat,distmat,hubmethod,Nreps)
% the probability that hubs remain after 100% rewire

n = size(mat,1) ;
% enfore und
mat = (mat + mat') ./ 2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% original hubs
[~,orighub] = get_hub_score_wei_und(mat,hubmethod) ;

rep_hubs = zeros(n,Nreps) ;
rep_jac = zeros(Nreps,1) ;
for ndx = 1:Nreps
    % disp(ndx)
    % full (100%) rewire using Gs
    [~,rewr_mat] = geombinsurr_partial(mat,distmat,1,20,'quantiles') ;
    % find the hubs
    [~,rep_hubs(:,ndx)] = get_hub_score_wei_und(rewr_mat,hubmethod) ;
    rep_jac(ndx) = jaccard(logical(orighub),logical(rep_hubs(:,ndx))) ;
end

% proportion of times hub is in same place
hubresil = mean(rep_hubs,2) ;
resiljac = mean(rep_jac) ; % overlap between orighubs and randomized hubs
resilvar = std(rep_jac) ; % variance of the overlap measure

%% calculate hubresil jaccard w/ orig @ different thresholds

hub_thresh = 0:0.01:1 ;
tmp_jacsim = zeros(length(hub_thresh),1) ;

for idx = 1:length(hub_thresh)  
    ttt = hub_thresh(idx) ;
    tmp_jacsim(idx) = jaccard(logical(hubresil>ttt),logical(orighub)) ;
end

hubresilarea = trapz(hub_thresh,tmp_jacsim) ;
