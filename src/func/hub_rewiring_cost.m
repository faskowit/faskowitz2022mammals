function [hubcostf, feedercostf, peripcostf,hubedged] = ...
    hub_rewiring_cost(mat,distmat,costmat,hubmethod,Nreps,prop_rand) 

if nargin < 4
   hubmethod = 'gollo2018' ; 
end

if nargin < 5
   Nreps = 1000 ;
end

if nargin < 6
   prop_rand = [ 0 ((5:5:100) ./ 100) ] ; % proprotion of randomizations 
end

if prop_rand(1) ~= 0 
   error('first value of prop_rand needs to be 0')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup stuff

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initalize stuff

hubcostf = zeros(length(prop_rand),2) ;
feedercostf = zeros(length(prop_rand),2) ;
peripcostf = zeros(length(prop_rand),2) ;

hubedged = zeros(length(prop_rand),1) ;

[~,orighubs] = get_hub_score_wei_und(mat,hubmethod) ;

% overallcost = sum(triunroll(costmat),'omitnan') ;

% initial hub cost
hubcostf(1,1) = sum(triunroll(costmat(orighubs,orighubs)),'omitnan') ;
feedercostf(1,1) = sum(costmat(~orighubs,orighubs),'all','omitnan') ;
peripcostf(1,1) = sum(triunroll(costmat(~orighubs,~orighubs)),'omitnan')  ;

% initalize hub cost reduction
hubedged(1) = density_und(mat(orighubs,orighubs)>0) ; 

%% 

for idx = 2:length(prop_rand)

    % disp(idx)
    currprop = prop_rand(idx) ;
    
    % calculate hub probability
    rep_hubcostfac = zeros(Nreps,1) ;
    rep_feedercostfac = zeros(Nreps,1) ;
    rep_periphcostfac = zeros(Nreps,1) ;
    rep_hubdens = zeros(Nreps,1) ;
    
    for ndx = 1:Nreps
        % disp(ndx) 
        % rewire using Gs
        [~,grewr_mat] = geombinsurr_partial(mat,distmat,currprop,50,'equalwidth') ;
        % find the hubs
        [~,hh] = get_hub_score_wei_und(grewr_mat,hubmethod) ;
        % get distances between new hubs
        
        rep_hubcostfac(ndx) = sum(triunroll(costmat(hh,hh)),'omitnan') ; 
        rep_feedercostfac(ndx) = sum(costmat(~hh,hh),'all','omitnan') ; 
        rep_periphcostfac(ndx) = sum(triunroll(costmat(~hh,~hh)),'omitnan') ; 
       
        rep_hubdens(ndx) = density_und(mat(hh,hh)>0) ;
        
    end
    
    hubcostf(idx,:) = [ mean(rep_hubcostfac) std(rep_hubcostfac) ] ;
    feedercostf(idx,:) = [ mean(rep_feedercostfac) std(rep_feedercostfac) ] ;
    peripcostf(idx,:) = [ mean(rep_periphcostfac) std(rep_periphcostfac) ] ;

    hubedged(idx) = mean(rep_hubdens) ;
end

