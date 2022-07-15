function [CN,R] = f_CN_corr_aprox(CIJ,nterms)

% Approximation of CN using Barnett et al. 
% FC is computed using Galan's linear model

% References
%   L. Barnett, C. L. Buckley and S. Bullock (2009)                      %
%   R. F. Galan, plos one (2008)

if nargin < 2
    nterms = 5;
end

N = size(CIJ,1);

% Analytical derivation of the covariance matrix
% parameters linear model

dt = 0.5;
alpha = 2;  
sigma = 1;
conn_factor = 1.0;
wij = 1;

W = CIJ.*wij;
W(isnan(W)) = 0;
normW = norm(W);
W = conn_factor*W/normW;
A = (1-alpha*dt)*eye(N) + W*dt;
A = (A+A')./2;

[VA, DA] = eig(A);
DA = diag(DA);
% spectral radius should be smaller than 1 -> stationarity
%ro_A = max(DA);  

Q = (sigma*dt)^2*eye(N); 
Qt = ((VA)^(-1))*Q*((VA')^(-1)); 

DA = DA(:,ones(N,1));
DA = 1 - (DA.*DA');

P = Qt./DA;
COV = VA*P*(VA');  

% correlation matrix
S = diag(1./sqrt(diag(COV)));
R = S*COV*S;
R = R-diag(diag(R));

m = 2;
mcn = zeros(1,nterms);

%[v d] = eig(R);
for term = 1:nterms
    mpowR = R^m;
    mcn(term) = ( ((-1)^m)*(m-1)/(m*(m+1)) )*trace(mpowR);
    m = m+1;
end

mcn = mcn*(N+1)/4;
CN = sum(mcn);

%err = trace(R^(m+1));




