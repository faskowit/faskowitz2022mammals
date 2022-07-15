function [hubs,hubprob,hubdarea,susceptind] = ...
    hub_sweep(mat,distmat,hubmethod,Nreps,proprand,altdist)
% from Gollo2018, get measures of hub fagility and measure the volatility
% of the hub distance embedding
%
% INPUTS:
%           mat:        connectivity data
%           distmat:    distance data
%           hubmethod:  method used to detect hubs, i.e. different hubscore
%                       definitions from previous papers. check out the
%                       function get_hub_score_... for more information
%           Nreps:      repetitions at each randomization proportion
%           proprand:   vector of randomization proportion values in the
%                       range 0-1; first value must be 0; ascending
%           altdist:    alternative distance mat for computation of
%                       inter-hub distance calculation.. otherwise, this
%                       dist just computed with distmat
%
% OUTPUTS:
%           hubs:       detected hubs on original data, according to hub
%                       method indicated
%           hubprob:    nodesXlen(proprand) matrix of proportion (prob)
%                       that a node is classified as a hub
%           hubdarea:   a measure of how the area of hubs changes over the 
%                       randomizations 
%           susceptind: the standard dev of the strengths at each level of 
%                       randomizations
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather inputs

if nargin < 3
   hubmethod = 'gollo2018' ; 
end

if nargin < 4
   Nreps = 100 ;
end

if nargin < 5
   proprand = [ 0 ((5:5:100) ./ 100) ] ; % proprotion of randomizations 
end

if nargin < 6
   altdist = distmat .* 1 ; 
end

if proprand(1) ~= 0 
   error('first value of prop_rand needs to be 0')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep vars

N = size(mat,1) ;
triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 

hubprob = zeros(N,length(proprand)) ; % init the hub rand mat
hubdarea = zeros(length(proprand),1) ; % init hub total distance
susceptind = zeros(length(proprand),1) ; % init susceptibility 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first level, no randomization

% first column is no randomization
% function [ hubscore, hubs] = get_hub_score_wei_und(mat,hubdeff,prcnthub)
[~,hubprob(:,1)] = get_hub_score_wei_und(mat,hubmethod) ;

% document the hubs, also for output
hubs = hubprob(:,1)==1 ;

% get initial hub distances
hubdarea(1) = sum(triunroll(altdist(hubs,hubs)),'omitnan') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do randomizations

timimg = 0 ; % debugging and stuff

for idx = 2:length(proprand)

    if timimg ; disp(idx) ; tic ; end %#ok<*UNRCH>
    
    currprop = proprand(idx) ;
    
    % calculate hub probability
    hubs_repout = zeros(N,Nreps) ;
    hdist_repout = zeros(Nreps,1) ;
    nstrengths = zeros(N,Nreps) ;
    for ndx = 1:Nreps
        % rewire using Gs
        [~,rewr_mat] = geombinsurr_partial(mat,distmat,currprop,20,'quantiles') ;
        % find the hubs
        [~,hubs_repout(:,ndx)] = get_hub_score_wei_und(rewr_mat,hubmethod) ;
        % get distances between new hubs
        hdist_repout(ndx) = sum(triunroll(altdist(~~hubs_repout(:,ndx),~~hubs_repout(:,ndx))),'omitnan') ;
        % get node strengths
        nstrengths(:,ndx) = sum(rewr_mat,2,'omitnan') ;
    end
    hubprob(:,idx) = mean(hubs_repout,2) ;
    hubdarea(idx) = mean(hdist_repout) ;
    susceptind(idx) = mean(std(nstrengths,[],2)) ;
    
    if timimg ; toc ; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hub volatility measures

% normalize against original total hub distance
hubdarea = hubdarea ./ hubdarea(1) ; 

% % get the area under the trapezoid to show hubs moving 'inward' 
% % small -> hubs compact quickly
% % large -> hubs distances don't change much
% hubdarea = trapz(proprand(:),hubdistnorm(:)) ;

% normalize suscept by max
susceptind = susceptind ./ max(susceptind) ;
