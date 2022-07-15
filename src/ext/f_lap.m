function [L,Ln,Lrw] = f_lap(adj)
% https://github.com/aiavenak/matlab_code/blob/master/f_lap.m

D = diag(sum(adj));
L = D - adj;

Ln = D^(-1/2)*L*D^(-1/2);

Lrw = D^(-1)*L;

end