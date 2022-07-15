function [meas] = get_comm_corr(wei,dist)

inv_wei = 1./wei ;
meas = nan(4,1) ;
m = logical(triu(ones(s+ize(wei,1)),1)) ;

% eff (inv shortest path)
try
    sp = 1./distance_wei_floyd(inv_wei) ; 
    meas(1) = corr(sp(m),dist(m),'type','s') ;  
catch
    warning('could not comput eff')
end

% -search information
try
    si = -1 .* search_information(wei,inv_wei,0) ;
    meas(2) = corr(si(m),dist(m),'type','s') ; 
catch
   warning('could not comput si') 
end

% diffusion eff (inv mfpt)
try 
    mt = 1./ mean_first_passage_time(wei) ; 
    mt = (mt + mt')./2 ; % make symmetric
    meas(3) = corr(mt(m),dist(m),'type','s') ; 
catch
   warning('could not compute diff eff') 
end

% comm
try
    co = communicability_wei(wei) ; 
    meas(4) = corr(co(m),dist(m),'type','s') ; 
catch
    warning('could not compute communicability')
end

