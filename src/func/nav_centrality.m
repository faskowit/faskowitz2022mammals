function [navcent] = nav_centrality(navdist,paths) 

n = size(navdist,1);
navcent = zeros(n,1) ;

for idx = 1:n
    for jdx = 1:n
        if isinf(navdist(idx,jdx)) || idx == jdx ; continue ; end     
        
        pa = paths{idx,jdx} ;
        navcent(pa(2:(end-1))) = navcent(pa(2:(end-1)))+1 ;
        
    end
end