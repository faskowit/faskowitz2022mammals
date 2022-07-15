function [hubs,hubfrag,hubprob,hubdista,susceptind] = ...
    hub_fragility_volatility(mat,distmat,hubmethod,Nreps,proprand,altdist)
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
%           hubfrag:    hub fragility; 1-> very fragile hubs, which
%                       dissappear quickly when perturbed
%           hubprob:    nodesXlen(proprand) matrix of proportion (prob)
%                       that a node is classified as a hub
%           hubdista:   a measure of how the area of hubs changes over the 
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
hubdista = zeros(length(proprand),1) ; % init hub total distance
susceptind = zeros(length(proprand),1) ; % init susceptibility 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first level, no randomization

% first column is no randomization
% function [ hubscore, hubs] = get_hub_score_wei_und(mat,hubdeff,prcnthub)
[~,hubprob(:,1)] = get_hub_score_wei_und(mat,hubmethod) ;

% document the hubs, also for output
hubs = hubprob(:,1)==1 ;

% get initial hub distances
hubdista(1) = sum(triunroll(altdist(hubs,hubs)),'omitnan') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do randomizations

timimg = 0 ; % debugging and stuff

for idx = 2:length(proprand)

    if timimg ; disp(idx) ; tic ; end %#ok<*UNRCH>
    
    currprop = proprand(idx) ;
    
    % calculate hub probability
    hubs_repout = zeros(N,Nreps) ;
    hd_repout = zeros(Nreps,1) ;
    nstrengths = zeros(N,Nreps) ;
    for ndx = 1:Nreps
        % rewire using Gs
        [~,rewr_mat] = geombinsurr_partial(mat,distmat,currprop,20,'quantiles') ;
        % find the hubs
        [~,hubs_repout(:,ndx)] = get_hub_score_wei_und(rewr_mat,hubmethod) ;
        % get distances between new hubs
        hd_repout(ndx) = sum(triunroll(altdist(~~hubs_repout(:,ndx),~~hubs_repout(:,ndx))),'omitnan') ;
        % get node strengths
        nstrengths(:,ndx) = sum(rewr_mat,2,'omitnan') ;
    end
    hubprob(:,idx) = mean(hubs_repout,2) ;
    hubdista(idx) = mean(hd_repout) ;
    susceptind(idx) = mean(std(nstrengths,[],2)) ;
    
    if timimg ; toc ; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hub volatility measures

% can do this outside of the function if you want to do adjustment
% % mean hub distances relative to the other nodes network
% relhub = mean(triunroll(distmat(hubs,hubs))) / mean(triunroll(distmat(~hubs,~hubs))) ;

% normalize against original total hub distance
hubdistnorm = hubdista ./ hubdista(1) ; 

% get the area under the trapezoid to show hubs moving 'inward' 
% small -> hubs compact quickly
% large -> hubs distances don't change much
hubdista = trapz(proprand(:),hubdistnorm(:)) ;

% normalize suscept by max
susceptind = susceptind ./ max(susceptind) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hub fragility meas

% now, get just the hubs, and look for when they pass 0.8
justhubs = hubprob(logical(hubs),:) ;
% and initalize the hub fragility index
hubfrag_tmp = zeros(size(justhubs,1),1) ;

% go through each node (row)
for idx = 1:size(justhubs,1)
    
    hubr = justhubs(idx,:) ; 
    hthr = hubr < 0.8 ; 
    crossind = find(hthr,1,'first') ; 
    if isempty(crossind) % if never goes below 0.8
        hubfrag_tmp(idx) = 0 ; % not fragile at all 
        continue 
    end   
    
    xcross = [ proprand(crossind-1) proprand(crossind) ] ;
    ycross = [ hubr(crossind-1) hubr(crossind) ] ;
    
    % slope (rise/run)
    m = (ycross(2) - ycross(1)) / (xcross(2) - xcross(1) ) ;
    % y intersect
    b = ycross(1) - (m*xcross(1)) ;
    % find point where 0.8 is crossed
    x0 = (0.8-b)/m; 
                        
    hubfrag_tmp(idx) = 1-x0 ;
                  
end

% output hub fragility as n length
hubfrag = nan(N,1) ;
hubfrag(logical(hubs)) = hubfrag_tmp ;
