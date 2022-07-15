function [ hubscore, hubs] = get_hub_score_bin_und(mat,hubdeff,prcnthub)
% various ways to determine which nodes might be hubs, but the binarized
% version of network functions, which makes some network functions redunant
% and makes some definitions slightly different than the original papers

if nargin<2
   hubdeff = 'betzel2014' ; 
end

if nargin<3
   prcnthub = .15 ; 
end

if (prcnthub <= 0 || prcnthub >= 1)
   error('prcnthub must be within (0,1)') 
end

n = size(mat,1) ;
% enfore und
mat = (mat+mat') ./ 2 ;
% and binarize
mat = weight_conversion(mat,'binarize') ;
lasthubind = upper(n*prcnthub) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(hubdeff)
    case 'betzel2014'
        % degrees
        d = degrees_und(mat) ;
%         % strength
%         s = strengths_und(mat) ;
        % efficiency
        e = efficiency_bin(mat,2) ; 
        % betweenness
        b = betweenness_bin(mat) ; 
        meas = [ d(:) e(:) b(:) ] ;
    case 'gollo2018'
        s = degrees_und(mat) ;
        meas = s(:) ;
    case 'perry2015'
        d = degrees_und(mat) ;
        b = betweenness_bin(mat) ; 
        sg = subgraph_centrality(mat) ;
        meas = [ d(:) b(:) sg(:) ] ;
    case 'swanson2016'
        d = degrees_und(mat) ;
%         s = strengths_und(mat) ;
        b = betweenness_bin(mat) ; 
        c = 1./sum(distance_bin(mat),2) ; % closeness cent
        meas = [ d(:) b(:) c(:) ] ;
    case 'sporns2007' 
        % this is a riff on paper, no formal hubscore given, the motif
        % finding takes a long time to run
        
        d = degrees_und(mat) ;
        b = betweenness_bin(mat) ; 
        c = 1./sum(distance_bin(mat),2) ; % closeness cent
        % frequency of motif 9
        [~,tmpf] = motif3struct_bin(mat) ;
        mf = tmpf(9,:) ; 
        meas = [ d(:) b(:) c(:) mf(:) ] ;
    case 'bassett2008' 
        % riff on the original method
        
        d = degrees_und(mat) ;
        b = betweenness_bin(mat) ; % do the inversion of weight
        c = 1./sum(distance_bin(mat),2) ; % closeness cent
        eg = eigenvector_centrality_und(mat) ;
        
        tmpm = zscore([ d(:) b(:) c(:) eg(:) ]) ;
        meas = tmpm .* double(tmpm>=2); % only count > 2 s.d. 
    case 'vandenheuvel2010' 
        % riff on original method
        
%         s = strengths_und(mat) ;
        d = degrees_und(mat) ; 
        b = betweenness_bin(mat) ; % do the inversion of weight
        c = 1./sum(distance_bin(mat),2) ; % closeness cent
        cc = -1 .* clustering_coef_bu(mat) ; % inver clustering coeff
        
        tmpm = [ d(:) b(:) c(:) cc(:) ] ;
        meas = tmpm .* double(tmpm > prctile(tmpm,67,1)) ; % top 33
    otherwise
        error('not a valid hubdeff option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rank_meas = abs(tiedrank(meas)-(n+1)) ; % make the score ascending (top=1st)
hubscore = mean(rank_meas,2) ; % average to get hubscore
[~,sort_hubscore] = sort(hubscore) ; % sort hub scor
hubs = false(n,1) ; 
hubs(sort_hubscore(1:lasthubind)) = 1 ; % output hubs

end
