function [C,Cthr] = get_cohesion_matrix(D, sym)
% https://www.pnas.org/content/119/4/e2003634119
% https://github.com/moorekatherine/pald/blob/main/R/pald_functions.R

if nargin < 2
    sym = 1 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(D,1) ;
C = zeros(size(D)) ;

for x = 1:(n-1)
    for y = ((x+1):n)
        
        dx = D(:,x) ;
        dy = D(:,y) ; 
        
        Uxy = find((dx <= D(y, x)) | (dy <= D(x, y))) ;
        
        wx = 1 * (dx(Uxy) < dy(Uxy)) + 0.5 * ( dx(Uxy) == dy(Uxy) ) ; 
        
        C(x, Uxy) = C(x, Uxy) + 1 ./ (length(Uxy)) .* wx' ; 
        C(y, Uxy) = C(y, Uxy) + 1 ./ (length(Uxy)) .* (1 - wx)' ;
        
    end
end

C = C ./(n-1) ;

if sym
   C = min(C,C') ; 
end

if nargout > 1
    thr = mean(diag(C)) / 2 ;
    Cthr = C .* 1 ;
    Cthr(C<thr) = 0 ;
end
