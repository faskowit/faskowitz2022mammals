function [rou,diffu] = scaled_routing_n_diffusion_und(mat,reps,rwiter,transform)

if nargin < 2
   reps = 10 ;
end

if nargin < 3
   rwiter = 1 ;  
end

if nargin < 4
   transform = 'inv' ; 
end

% clean up mat
gg = get_components(mat) ;
ggmode = mode(gg) ;
if any(gg~=ggmode(1))
   warning('more than one component, will fix') 
   mat = mat(gg==ggmode(1),gg==ggmode(1)) ;
end

N = size(mat,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rou_emp = rout_efficiency(mat,transform);
diffu_emp = diffusion_efficiency(mat) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latmio_nets = zeros(N,N,reps) ;
randmio_nets = zeros(N,N,reps) ;

% generate 
for idx = 1:reps
    % disp(idx)
    
    latmio_nets(:,:,idx) = latmio_und_connected(mat,rwiter) ;
    randmio_nets(:,:,idx) = randmio_und_connected(mat,rwiter) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rou_lat = arrayfun(@(ind) rout_efficiency(latmio_nets(:,:,ind),transform),1:reps) ;
diffu_lat = arrayfun(@(ind) diffusion_efficiency(latmio_nets(:,:,ind)),1:reps) ;

rou_rand = arrayfun(@(ind) rout_efficiency(randmio_nets(:,:,ind),transform),1:reps) ;
diffu_rand = arrayfun(@(ind) diffusion_efficiency(randmio_nets(:,:,ind)),1:reps) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rou = (rou_emp - mean(rou_lat)) / ( mean(rou_rand) - mean(rou_lat) ) ;
diffu = (diffu_emp - mean(diffu_lat)) / ( mean(diffu_rand) - mean(diffu_lat) ) ;
