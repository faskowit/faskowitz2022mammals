function Dab = lap_spect_dist(aij,bij,nComps,minkp)
% https://github.com/aiavenak/matlab_code/blob/master/f_lap_dist.m
% calculate the edit distance between two graphs A and B
% reference: Thomas Thorne and Michael P. H. Stumpf
%            evolution Graph spectral analysis of 
%            protein interaction network
%            J. R. Soc. Interface 2 May 2012

if nargin < 3
   nComps = 0 ; % if 0, all components
end

if nargin < 4
   minkp = 2 ; % 2 -> euclidean distance 
end

% [d1,d2] = size(aij);
% [b1,b2] = size(bij);
% if (d1 == b1) && (d2 == b2)
%     if d1 == d2  % aij, bij are matrices
        % get laplacian and eigenvalues
        [~,Lna] = f_lap(aij);
        [~,Lnb] = f_lap(bij);
        
        if nComps > 0
            e_opts.tol = 1e-8;
            e_opts.disp = 0;
            e_opts.isreal = 1;
            
            [~,Va] = eigs(double(Lna),nComps+1,'smallestreal',e_opts);
            [~,Vb] = eigs(double(Lnb),nComps+1,'smallestreal',e_opts);      
        else
            [~,Va] = eig(Lna);
            [~,Vb] = eig(Lnb);
        end
            
        Va = diag(Va);
        Vb = diag(Vb);

        Va = sort(Va,'ascend');
        Vb = sort(Vb,'ascend');

%     else  % aij, bij are vectors of eigenvalues
%         Va = aij;
%         Vb = bij;
%     end
    
    Dab = pdist([Va' ; Vb'],'minkowski',minkp) ;
    
% else
%     warning('sizes do not match');
%     Dab = 0/0;
% end

end
