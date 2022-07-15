function [qcontributions,q,B] = fcn_get_q(W,M,B,gamma)
%%FCN_GET_Q         get modularity of partition
%
%   Q = FCN_GET_Q(CI,AIJ) given a partition CI and adjacency matrix AIJ,
%       compute the modularity score of the partition. This script is
%       applicable to only symmetric networks.
%
%   Inputs:     W,              connectivity matrix (dimensions [n x n])
%               M,              community labels (dimensions [n x 1])
%               B,              modularity type
%               GAMMA,          resolution parameter
%
%   Outputs:    Qcontributions, contributions to total Q by each community
%                               (dimensions [N_c x 1])
%               Q,              total modularity of the partition (1 x 1)
%               B,              modularity matrix used to calculate
%
%   Example usage:
%
%   
%   load ts;                                    % load time series data
%   load system_labels;                         % load system assignments
%   W = corr(ts);                               % compute correlation matrix
%   Qc = fcn_get_q(W,lab,'negative_asym',1);    % get contributions from each system using negative, asymmetric modularity (Rubinov & Sporns, 2011)
%
%   Richard Betzel, Indiana University, 2012

W=double(W);                                % convert to double format
s=sum(sum(W));                              % get sum of edges

if ~exist('B','var') || isempty(B)
    type_B = 'modularity';
elseif ischar(B)
    type_B = B;
else
    type_B = 0;
    if exist('gamma','var') && ~isempty(gamma)
        warning('Value of gamma is ignored in generalized mode.')
    end
end
if ~exist('gamma','var') || isempty(gamma)
    gamma = 1;
end

if strcmp(type_B,'negative_sym') || strcmp(type_B,'negative_asym')
    W0 = W.*(W>0);                          %positive weights matrix
    s0 = sum(sum(W0));                      %weight of positive links
    B0 = W0-gamma*(sum(W0,2)*sum(W0,1))/s0; %positive modularity
    
    W1 =-W.*(W<0);                          %negative weights matrix
    s1 = sum(sum(W1));                      %weight of negative links
    if s1                                   %negative modularity
        B1 = W1-gamma*(sum(W1,2)*sum(W1,1))/s1;
    else
        B1 = 0;
    end
elseif min(min(W))<-1e-10
    err_string = [
        'The input connection matrix contains negative weights.\nSpecify ' ...
        '''negative_sym'' or ''negative_asym'' objective-function types.'];
    error(sprintf(err_string))              
end
if strcmp(type_B,'potts') && any(any(W ~= logical(W)))
    error('Potts-model Hamiltonian requires a binary W.')
end

if type_B
    switch type_B
        case 'modularity';      B = (W-gamma*(sum(W,2)*sum(W,1))/s)/s;
        case 'potts';           B =  W-gamma*(~W);
        case 'negative_sym';    B = B0/(s0+s1) - B1/(s0+s1);
        case 'negative_asym';   B = B0/s0      - B1/(s0+s1);
        otherwise;              error('Unknown objective function.');
    end
else                            % custom objective function matrix as input
    B = double(B);
    if ~isequal(size(W),size(B))
        error('W and B must have the same size.')
    end
end

h = dummyvar(M);
qcontributions = diag(h'*(B*h));
q = sum(qcontributions);