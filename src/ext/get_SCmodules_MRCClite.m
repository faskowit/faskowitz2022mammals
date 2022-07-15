function [ciu,Aall,anull,A,ciall] = get_SCmodules_MRCClite(K,sam1,sam2,maxC,r1,r2,tau)
% 
%   INPUTS:
%       K:      network data
%       sam1:   samples communities first loop
%       sam2:   samples communities second loop
%       maxC:   max num communities
%       r1:     initial gamma low
%       r2:     initial gamma high
%       tau:    tau parameter for consensus modules at the end 
%
%   OUTPUTS:
%       ciu:    consensus community partition
%       Aall:   agreement matrix
%       anull:  analytic null
%       A:      agree - null
%       ciall:  community assignments

if nargin < 7
   tau = 0.0 ; 
end

% performs a 'lite' version of MRCC, obtaining a CC (Aall) matrix and a set of
% consensus modules

N = size(K,1);

% FIRST PASS - find gamma range where modules range from 2 to (just below) maxC;
% initial gamma range
gam_range = logspace(r1,r2,sam1);
G = length(gam_range);

ciall = zeros(N,G);
for g=1:length(gam_range)
    if ~mod(g,100) ; disp([ num2str(g) ' of ' num2str(length(gam_range))]) ; end
    [ci, ~] = community_louvain(K,gam_range(g),[]);
    ciall(:,g) = ci;
end
% identify the limits of the range
ff = find((max(ciall)>1)&(max(ciall)<maxC));
g1 = log10(gam_range(ff(1)));
g2 = log10(gam_range(ff(end)));
disp([num2str(g1),' ',num2str(g2)])

% SECOND PASS - use the gamma range determined in first pass
% collect all partitions within that range
gam_range = logspace(g1,g2,sam2);
G = length(gam_range);

ciall = zeros(N,G);
for g=1:G
    if ~mod(g,1000) ; disp([ num2str(g) ' of ' num2str(G)]) ; end
    [ci,~] = community_louvain(K,gam_range(g),[]);
    ciall(:,g) = ci;
end

% exclude spurious partitions that are outside of the range 2 to maxC
numm = max(ciall);
use = find((numm>1)&(numm<maxC));
ciall = ciall(:,use);

% co-classification matrix
Aall = agreement(ciall)./length(use);

% analytic null
anull = 0;
for cnum = 1:length(use)
    anull = anull+sum((histcounts(ciall(:,cnum))./N).*((histcounts(ciall(:,cnum))-1)./(N-1)));
end
anull = anull/length(use);
A = Aall - anull;

% consensus clustering
ciu = consensus_und(A,tau,100);
