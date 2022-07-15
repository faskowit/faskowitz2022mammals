function [hubfrag80,hubfragAUC] = get_hub_fragility(hubs,hubprob,proprand) 
% 
% get the hub fragility, given hubs and hub probabilities. use function in
% tandem with hub_sweep
%
% INPUTS:
%   hubs:       Nx1 vector of hubs
%   hubprob:    Nxlength(proprand) matrix of hub probabilities
%   proprand:   1xN vector of randomization proportions used when
%               calculating the hubprob data
% OUTPUTS:
%   hubfrag80:    fragility of original hubs, based on crossing .8 prob.
%   hubfragAUC:   area under curve across all randomizations
%

if length(hubs) ~= size(hubprob,1)
   error('make sure hubs vec same node dimension as hubprob') 
end

if size(hubprob,2) ~= length(proprand)
   error('make sure hubprob has same 2nd dimension as proprand') 
end

if proprand(1) ~= 0 
   error('first value of proprand should be 0')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(hubs) ;
fragility_thr = 0.8 ; % given by paper

% get just the hubs, and look for when they pass 0.8
justhubs = hubprob(logical(hubs),:) ;
% and initalize the hub fragility index
hubfrag80_tmp = zeros(size(justhubs,1),1) ;
hubfragAUC_tmp = nan(size(justhubs,1),1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go through each hub node (row)
for idx = 1:size(justhubs,1)
    
    hubr = justhubs(idx,:) ; 
    
    % fragility auc
    hubfragAUC_tmp(idx) = 1-trapz(proprand,hubr) ;
    
    % fragility 80 
    hthr = hubr < fragility_thr ; 
    crossind = find(hthr,1,'first') ; 
    if isempty(crossind) % if never goes below 0.8
        hubfrag80_tmp(idx) = 0 ; % not fragile at all 
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
                        
    hubfrag80_tmp(idx) = 1-x0 ;
                  
end

% output hub fragility as n length
hubfrag80 = nan(N,1) ;
hubfrag80(logical(hubs)) = hubfrag80_tmp ;

hubfragAUC = nan(N,1) ;
hubfragAUC(logical(hubs)) = hubfragAUC_tmp ;

