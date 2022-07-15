function [ out_struct ] = readin_mammal_dat(filename) 

if isfile(filename)
    ll = load(filename) ;
else
    error('file does not exist') 
end

% spatial
out_struct.dist_mat = squareform(pdist(ll.cents)) ; 
out_struct.leng_mat = ll.LmatT ;
out_struct.node_adj = get_seg_spatial_adj(ll.seg,1:200,6) ;
out_struct.cents = ll.cents ;

% vol 
out_struct.segs_vox = arrayfun(@(x)sum(ll.seg(:) == x), 1:200)' ;
geo_means = (out_struct.segs_vox*out_struct.segs_vox').^0.5 ; 

% connectivity
out_struct.con_mat = ll.conmatT ;
out_struct.con_mat_gn = ll.conmatT ./ geo_means ;
