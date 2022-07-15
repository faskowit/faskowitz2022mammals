function [adj_mat] = get_seg_spatial_adj(seg,vals,conn)
% given 3-D segmentation, find regions that are adjacent
% uses the function findNeighbours found here: 
% https://www.mathworks.com/matlabcentral/fileexchange/68549-findneighbours

if nargin < 3
   conn = 26 ; 
end

if ~ismember(conn,[ 6 18 26 ])
    error('conn needs to be either 6 18 26')
end

n_vals = length(vals) ;
uniq_in_seg = unique(nonzeros(seg)) ;
adj_mat = zeros(n_vals,length(uniq_in_seg)) ;
seg_size = size(seg)  ;

for idx = 1:n_vals
    
   % disp(idx) 
   
   curr_val = vals(idx) ;
   
   % get inds where seg is curr_val
   curr_locs = find(seg==curr_val) ; 
   
   curr_nei = [] ; 
   for jdx = 1:length(curr_locs)
      
       nn = findNeighbours(curr_locs(jdx),seg_size,conn) ; 
       curr_nei = [ curr_nei ; nonzeros(seg(nn)) ] ;
   end
   
   curr_nei(curr_nei==curr_val) = [] ;
   uniq_nei = unique(curr_nei) ;
   
   adj_vals = arrayfun(@(x)sum(curr_nei == x), uniq_nei) ;
   
   % map to ind, so as not to assume the vals are the inds
   [uniq_ind,~] = ismember(uniq_in_seg,uniq_nei);
   
   adj_mat(curr_val,uniq_ind) = adj_vals ;
   
end



