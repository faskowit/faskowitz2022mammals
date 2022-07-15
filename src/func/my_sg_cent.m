function [sgc] = my_sg_cent(W)
%inputs
%           W      weighted connection matrix
%
%outputs
%           co      subgraph centrality
%=================================================

B = sum(W,2);
C = diag(B);
D = C^(-(1/2));
E = D * W * D;
F = expm(E);
sgc = diag(F) ;

end