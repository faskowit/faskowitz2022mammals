function [net_dist] = netpd_rewire_dist(mat,distmat,rewiremeth,Nreps,prop_rand,nbins)
% inspired by Gollo, measure network overall distance using the network
% portait divergence

if nargin < 3 
    rewiremeth = 'randmio' ; % func from BCT
end

if nargin < 4
   Nreps = 100 ;
end

if nargin < 5
   prop_rand = [ 0 ((5:5:100) ./ 100) ] ; % proprotion of randomizations 
end

if nargin < 6
   nbins = 50 ; 
end

if prop_rand(1) ~= 0 
   error('first value of prop_rand needs to be 0')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup rewire func

switch lower(rewiremeth)
    case {'geombinsurr_gs','gs'} 
        rewirefunc = @(cp_) geombinsurr_gs(mat,distmat,cp_,20,'quantiles') ;
    case {'geombinsurr_gss','gss'} 
        rewirefunc = @(cp_) geombinsurr_gss(mat,distmat,cp_,20,'quantiles') ;
    case {'geombinsurr_gw','gw'} 
        rewirefunc = @(cp_) geombinsurr_partial(mat,distmat,cp_,20,'quantiles') ;
    case {'shuffwei','sw'}
        rewirefunc = @(cp_) shuffwei_partial(mat,cp_) ;
    case {'randmio','rm'}
        rewirefunc = @(cp_) randmio_und_connected(mat,cp_) ;
    otherwise
        error('Error: rewire method not valid: %s',rewiremeth)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initalize stuff

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 

% get initial portrait 
D1 = distance_wei_floyd(1./mat) ;
min_D1 = min(triunroll(D1)) ; 
max_D1 = max_noinf(triunroll(D1)) ; % noinf to account for unreachable nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do geometric randomizations

net_dist = zeros(length(prop_rand),2) ;

for idx = 2:length(prop_rand)

%     disp(idx)
    currprop = prop_rand(idx) ;
    
    % calculate rewiring distanc
    reps_d = zeros(Nreps,1) ;
    
    for ndx = 1:Nreps % ; disp(ndx) 
        
        % rewire 
        rewr_net = rewirefunc(currprop) ;

        % get distance mats from rewired matrices 
        D2 = distance_wei_floyd(1./rewr_net) ;

        % get edge bins
        eb=linspace(min([min_D1 min(triunroll(D2))]), ...
            max([ max_D1 max_noinf(triunroll(D2)) ] )+eps,nbins+1);

        % portraits, with common bins
        B1 = netpd_portrait_wei(D1,eb,'alreadydistance') ;
        B2 = netpd_portrait_wei(D2,eb,'alreadydistance') ;

        % the network portait divergence
        reps_d(ndx) = netpd_divergence(B1,B2) ;
    end
    
    net_dist(idx,:) = [ mean(reps_d) std(reps_d) ] ;
end

end

function m = max_noinf(dat)
% get max with no inf
    m = max(dat(~isinf(dat)));
end

function [rewnet] =  geombinsurr_gs(varargin)
% need function to capture specific output
    [~,rewnet] = geombinsurr_partial(varargin{:}) ;
end

function [rewnet] =  geombinsurr_gss(varargin)
% need function to capture specific output
    [~,~,rewnet] = geombinsurr_partial(varargin{:}) ;
end

