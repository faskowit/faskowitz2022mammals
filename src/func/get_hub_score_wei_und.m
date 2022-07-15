function [ hubscore, hubs, hcorr ] = get_hub_score_wei_und(mat,hubdeff,hubtype,prcnthub)
% various ways to determine which nodes might be hubs

if nargin<2
   hubdeff = 'betzel2014' ; 
end

if nargin<3
    hubtype = 'rank' ;
end

if strcmp(hubtype,'rank')
    hubrank = 1 ;
elseif strcmp(hubtype,'zscore')
    hubrank = 0 ;
else
    error('hubtype must be "rank" or "zscore"') 
end

if nargin<4
   prcnthub = .15 ; 
end

if (prcnthub <= 0 || prcnthub >= 1)
   error('prcnthub must be within (0,1)') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(mat,1) ;
% enfore und
mat = (mat + mat') ./ 2 ;
lasthubind = ceil(n*prcnthub) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(hubdeff)
    case 'mydef'
        % degrees
        d = degrees_und(mat) ;
        % strength
        s = strengths_und(mat) ;
        % efficiency
        cc = clustering_coef_wu(mat) ; 
        % betweenness
        b = my_btwn_cent(mat,'inv') ; 
        % closeness
        c = my_close_cent(mat) ; 
        meas = [ d(:) s(:) cc(:) b(:) c(:) ] ;
    case 'mydef2'
        % strength
        s = strengths_und(mat) ;
        % degree 
        d = degrees_und(mat) ;
        % subgraph centrality
        sg = my_sg_cent(mat) ;  
        % between
        b = my_btwn_cent(mat) ; 
        meas = [ s(:) d(:) sg(:) b(:) ] ;
    case 'betzel2014'
        % degrees
        d = degrees_und(mat) ;
        % strength
        s = strengths_und(mat) ;
        % efficiency
        e = my_efficiency_wei(mat,2,'inv') ; 
        % betweenness
        b = my_btwn_cent(mat,'inv') ; 
        meas = [ d(:) s(:) e(:) b(:) ] ;
    case 'gollo2018'
        s = strengths_und(mat) ;
        meas = s(:) ;
    case 'perry2015'
        d = degrees_und(mat) ;
        b = my_btwn_cent(mat) ; 
        sg = my_sg_cent(mat) ;
        meas = [ d(:) b(:) sg(:) ] ;
    case 'swanson2016'
        d = degrees_und(mat) ;
        s = strengths_und(mat) ;
        b = my_btwn_cent(mat) ;
        c = my_close_cent(mat) ; 
        meas = [ d(:) s(:) b(:) c(:) ] ;
    case 'sporns2007' 
        % this is a riff on paper, no formal hubscore given, the motif
        % finding function takes a long time to run
        
        d = degrees_und(mat) ;
        b = my_btwn_cent(mat) ; 
        c = my_close_cent(mat) ; 
        % frequency of motif 9
        [~,~,tmpf] = motif3struct_wei(weight_conversion(mat,'normalize')) ;
        mf = tmpf(9,:) ; 
        meas = [ d(:) b(:) c(:) mf(:) ] ;
    case 'bassett2008' 
        % riff on the original method
        
        d = degrees_und(mat) ;
        b = my_btwn_cent(mat) ; 
        c = my_close_cent(mat,'inv') ; 
        eg = eigenvector_centrality_und(mat) ;
        
        tmpm = zscore([ d(:) b(:) c(:) eg(:) ]) ;
        meas = tmpm .* double(tmpm>=2); % only count > 2 s.d. 
    case 'vandenheuvel2010' 
        % riff on original method
        
        s = strengths_und(mat) ;
        b = my_btwn_cent(mat) ; 
        c = my_close_cent(mat,'inv') ; 
        cc = -1 .* clustering_coef_wu(mat) ; % invert clustering coeff
        
        tmpm = [ s(:) b(:) c(:) cc(:) ] ;
        meas = tmpm .* double(tmpm > prctile(tmpm,67,1)) ; % top 33
    otherwise
        error('not a valid hubdeff option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hubrank
    % hubscore based on rankings 
    rank_meas = abs(tiedrank(meas)-(n+1)) ; % make the score ascending (top=1st)
    hubscore = mean(rank_meas,2) ; % average to get hubscore
    [~,sort_hubscore] = sort(hubscore) ; % sort hub score
    hubs = false(n,1) ; 
    hubs(sort_hubscore(1:lasthubind)) = 1 ; % output hubs
else
    % hubscore based on zscore
    hubscore = mean(zscore(meas),2) ;
    hubs = hubscore>=(norminv(1-prcnthub)) ; % output hubs
end

% correlation between hubscore measures
hcorr = corr(meas,'type','s') ; % spearman correlation btwn hub meas

end % end of main function

% add the better efficienty function
function E = my_efficiency_wei(W, local, transdist)
% slightly modified copy of efficiency function from BCT; added transdist
% and added distance floyd to speed things up drastically
% J.F 

%Modification history
% 2011: Original (based on efficiency.m and distance_wei.m)
% 2013: Local efficiency generalized to directed networks
% 2017: Added the modified local efficiency and updated documentation.

if nargin < 3
    transdist = 'inv' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(W);                                              % number of nodes
ot = 1 / 3;                                                 % one third

L = W .* 1; % make copy                                     % connection-length matrix
A = W > 0;                                                  % adjacency matrix

switch transdist
    case 'inv'
        L(A) = 1 ./ L(A);
    case 'log'
        L(A) = -log(L(A)./ (max(L(A))+eps)) ; % +eps makes no 0 distances
    otherwise
        error('not a valid transdist') 
end

A = double(A);

if exist('local','var') && local                            % local efficiency
    E = zeros(n, 1);
    cbrt_W = W.^ot;
    switch local
        case 1
            for u = 1:n
                V  = find(A(u, :) | A(:, u).');             % neighbors
                sw = cbrt_W(u, V) + cbrt_W(V, u).';       	% symmetrized weights vector
                di = distance_inv_wei(L(V, V));             % inverse distance matrix
                se = di.^ot + di.'.^ot;                     % symmetrized inverse distance matrix
                numer = (sum(sum((sw.' * sw) .* se)))/2;   	% numerator
                if numer~=0
                    sa = A(u, V) + A(V, u).';              	% symmetrized adjacency vector
                    denom = sum(sa).^2 - sum(sa.^2);        % denominator
                    E(u) = numer / denom;                   % local efficiency
                end
            end
        case 2
            cbrt_L = L.^ot;
            for u = 1:n
                V  = find(A(u, :) | A(:, u).');            	% neighbors
                sw = cbrt_W(u, V) + cbrt_W(V, u).';       	% symmetrized weights vector
                di = distance_inv_wei(cbrt_L(V, V));      	% inverse distance matrix
                se = di + di.';                             % symmetrized inverse distance matrix
                numer=(sum(sum((sw.' * sw) .* se)))/2;      % numerator
                if numer~=0
                    sa = A(u, V) + A(V, u).';             	% symmetrized adjacency vector
                    denom = sum(sa).^2 - sum(sa.^2);        % denominator
                    E(u) = numer / denom;                 	% local efficiency
                end
            end
    end
else
    di = distance_inv_wei(L);
    E = sum(di(:)) ./ (n^2 - n);                         	% global efficiency
end

function D=distance_inv_wei(L)
% inverted distance using faster floyd
    
D = distance_wei_floyd(L) ;
D=1./D;                                                     % invert distance
D(1:(size(L,1)+1):end)=0;

end

end % end of my efficiency

% added better betweenness centrality
function btwncent = my_btwn_cent(W,transdist)
% using floyd. do the inverse weight inside func

if nargin < 2
    transdist = 'inv' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch transdist
    case 'inv'
        L = 1 ./ W;
    case 'log'
        L = -log(W./ (max(W)+eps)) ; % +eps makes no 0 distances
    otherwise
        error('not a valid transdist') 
end

n = size(W,1) ;
[~,hops_,paths_] = distance_wei_floyd(L) ;
btwncent = zeros(n,1) ; 

for idx = 1:n
    for jdx = 1:n
        if idx >= jdx ; continue ; end
        
        pa = retrieve_shortest_path(idx,jdx,hops_,paths_) ;
        btwncent(pa(2:(end-1))) = btwncent(pa(2:(end-1)))+2 ;
    end
end

end % end of my betweenness

function cc = my_close_cent(W,transdist)
% using floyd. do the inverse weight inside func

if nargin < 2
    transdist = 'inv' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch transdist
    case 'inv'
        L = 1 ./ W;
    case 'log'
        L = -log(W./ (max(W)+eps)) ; % +eps makes no 0 distances
    otherwise
        error('not a valid transdist') 
end

d = distance_wei_floyd(L) ;
d(isinf(d)) = nan ; % if there are infs... make into nan
dd = sum(d,2,'omitnan') ;
dd(dd==0) = eps ; % if there are 0's (row of nan)... make tiny tiny number
cc = 1./dd(:) ;

end % end of closeness

function sgc = my_sg_cent(W)
%inputs
%           W      weighted connection matrix
%
%outputs
%           co      subgraph centrality
%=================================================

B = sum(W,2);
C = diag(B);
D = C^(-(1/2));
E = D * W * D;
F = expm(E);
sgc = diag(F) ;

end

