clear all
close all
clc

load('../../data/Streamline_60_group.mat');
N = 40;
D = 0.25;
M = round(D*N*(N - 1)/2);
w = normrnd(0.5,0.15,M,1);
w(w < 0) = -w(w < 0);

K = 4;
G = N/K;
p = 0.8;
A = zeros(N);
ci = repmat(1:K,[G,1]);
ci = ci(:);

for i = 1:M
    [u,v] = find(triu(A == 0,1));
    idx = (v - 1)*N + u;
end

% 
% A = W.*G;
% m = 11;
% maxSamp = round(logspace(0,2,m));
% nreps = 10;
% 
% cn = zeros(nreps,m);
% for i = 1:m
%     for irep = 1:nreps
%         [cn(irep,i),err,R] = f_CN_smp(A,maxSamp(i));
%         imagesc(cn); drawnow;
%     end
% end