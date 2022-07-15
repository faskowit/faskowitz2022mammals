function [Gw,Gs,Gss]=geombinsurr_partial(W,D,frac,nbins,bintype)
% [Gw,Gs,Gss]=geombinsurr_partial(W,D,frac,nbins,[bintype])
%
% Surrogate random graphs that approximately preserve the relationship
% between weight and distance by first binning the edges by distance and
% then shuffling within the bins. Randomizes only a fraction frac of the
% edges.
%
% Inputs
% W         = weighted connectivity matrix
% D         = physical distance between nodes
% frac      = scalar 0<frac<=1, shuffle this fraction of edges
% nbins     = number of bins
% bintype   = 'quantiles' for bins containing equal numbers of edges
%             or 'equalwidth' for bins of equal width
%             (optional, default='equalwidth')
%
% Outputs
% Gw       = random surrogate preserving weights but not strengths
% Gs       = preserves strength distribution but not sequence
% Gss      = preserves strength sequence
%
% Note: Gs and Gss are generated assuming that W is undirected.
%
% References: Gollo et al. (2018) Nature Neuroscience.
%             Roberts et al. (2016) NeuroImage 124:379-393.
%
% J Roberts, M Breakspear, L Gollo
% QIMR Berghofer, 2018

if nargin<5
    bintype='equalwidth';
end

%Preliminaries
N = size(W,1);
Gw = zeros(N,N); 

%Check if directed
drct=1;
if max(max(W-W'))==0, drct=0; end

W(1:N+1:end)=0;         %Ensure there's no self connections
if drct==0              %Only do one triangle if undirected, then replicate
    W = triu(W);
end

nz= find(W);            %Nonzero entries
w = W(nz);              %Vector of non-zero connections
d = D(nz);              %Corresponding distances

%Set up bins by distance
switch bintype
    case 'quantiles'
        %bin by quantiles (bins contain equal numbers of edges)
        if nbins>2
            binedges=[min(d) quantile(d,nbins-1) max(d)+eps]; 
        elseif nbins==2
            binedges=[min(d) median(d) max(d)+eps];
        else
            binedges=[min(d) max(d)+eps];
        end
    case 'equalwidth'
        %bin into fixed-width bins 
        binedges=linspace(min(d),max(d)+eps,nbins+1);
    otherwise
        error('invalid bintype: %s (must be ''quantiles'' or ''equalwidth'')',bintype)
end

nedges=length(nz);
subsetsize=round(frac*nedges);
shuffsubseti=randsample(nedges,subsetsize); % choose a random subset
dsubset=d(shuffsubseti);
wsubset=w(shuffsubseti);

% bin the edges
[counts,inds]=histc(dsubset,binedges); % bin just the random subset of edges
for j=1:nbins
    % find which edges are in a given bin
    whichw=find(inds==j);  
    if counts(j)>0
        % shuffle all edges within the bin
        wsubset(whichw)=wsubset(whichw(randperm(counts(j))));  
    end
end

w(shuffsubseti)=wsubset; % put shuffled subset back into full set of weights

Gw(nz)=w;

if drct==0
    Gw=(Gw+Gw');
    W=W+W';
end

if nargout > 1 % Gs output
    %Adjust node strengths - NOTE: assumes W is undirected
    strengthsW=sum(W);     % original node strengths
    strengthsGw=sum(Gw); % new node strengths
    %Re-order the old strengths to match new random order in Gw
    strengthsGs=rankreorder(strengthsW,strengthsGw);
    %Adjust strengths to give both original and random strength sequences
    Gs=strengthcorrect(Gw,strengthsGs); % orig strengths in new sequence
    if nargout > 2 % Gss output 
        Gss=strengthcorrect(Gw,strengthsW);  % orig strengths in orig sequence
    end
end

end

function out=rankreorder(x,scaffold)
% reorder vector x to have same rank order as vector scaffold
y(:,1)=scaffold;
y(:,2)=1:length(x);
y=sortrows(y,1);
y(:,1)=sort(x);
y=sortrows(y,2);
out=y(:,1);
end
