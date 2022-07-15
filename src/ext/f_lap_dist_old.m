function Dab = f_lap_dist(aij,bij,nComps)
% https://github.com/aiavenak/matlab_code/blob/master/f_lap_dist.m
% calculate the edit distance between two graphs A and B
% reference: Thomas Thorne and Michael P. H. Stumpf
%            evolution Graph spectral analysis of 
%            protein interaction network
%            J. R. Soc. Interface 2 May 2012

if nargin < 3
   nComps = 0 ; 
end

[d1,d2] = size(aij);
[b1,b2] = size(bij);
if (d1 == b1) && (d2 == b2)
    if d1 == d2  % aij, bij are matrices
        % get laplacian and eigenvalues
        [La,Lna] = f_lap(aij);
        [Lb,Lnb] = f_lap(bij);
        
        if nComps > 0
            % href="https://brainspace.readthedocs.io/en/latest/pages/matlab_doc/support_functions/laplacian_eigenmaps.html">ReadTheDocs</a>.

            % Construct eigenmaps (solve Ly = lambda*Dy)
            % disp('Constructing Eigenmaps...');
            tol = 1e-10;
            options.disp = 0;
            options.isreal = 1;

            [~, Va] = eigs(double(La), double(diag(sum(aij))), nComps + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
            Va = diag(Va);
            Va = sort(Va, 'ascend');
            Va = Va(2:nComps + 1);

            [~, Vb] = eigs(double(Lb), double(diag(sum(bij))), nComps + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
            Vb = diag(Vb);
            Vb = sort(Vb, 'ascend');
            Vb = Vb(2:nComps + 1);
             
        else
            [~,Va] = eig(Lna);
            [~,Vb] = eig(Lnb);
            
            Va = diag(Va);
            Vb = diag(Vb);
            
            Va = sort(Va,'ascend');
            Vb = sort(Vb,'ascend');
        end

    else  % aij, bij are vectors of eigenvalues
        Va = aij;
        Vb = bij;
    end
    
    Dab = sum((Va-Vb).^2);
    
else
    warning('sizes do not match');
    Dab = 0/0;
end

end
