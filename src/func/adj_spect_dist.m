function Dab = adj_spect_dist(aij,bij,nComps,minkp)
% adapted from:
% https://github.com/aiavenak/matlab_code/blob/master/f_lap_dist.m

if nargin < 3
   nComps = 0 ; % if 0, all components
end

if nargin < 4
   minkp = 2 ; % 2 -> euclidean distance 
end

if nComps > 0
    e_opts.tol = 1e-8;
    e_opts.disp = 0;
    e_opts.isreal = 1;

    [~,Va] = eigs(double(aij),nComps,'largestreal',e_opts);
    [~,Vb] = eigs(double(bij),nComps,'largestreal',e_opts);      
else
    [~,Va] = eig(aij);
    [~,Vb] = eig(bij);
end

Va = diag(Va);
Vb = diag(Vb);

Va = sort(Va,'descend');
Vb = sort(Vb,'descend');

Dab = pdist([Va' ; Vb'],'minkowski',minkp) ;

end
