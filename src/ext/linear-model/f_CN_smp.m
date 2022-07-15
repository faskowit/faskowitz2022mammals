function [cn,err,R] = f_CN_smp(CIJ,maxSamp)

% Approximation of CN using "sampling method"
% FC is computed using Galan's linear model

% References
%   L. Barnett, C. L. Buckley and S. Bullock (2009)                      %
%   R. F. Galan, plos one (2008)

if nargin < 2
    maxSamp = 500;
end

N = size(CIJ,1);

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

err = 0;
cn = 0;
n = N;

logd = logdet(COV);
if isinf(logd)
	fprintf(2,'ERROR (CN): bad covariance matrix\n');
	err = 2;
	return;
end

% suppress warning resulting from large nsamp value 
% when nchoosek is too big
warning off MATLAB:nchoosek:LargeCoefficient

for k = 1:n-1
    
    nsamp = ceil(0.1*nchoosek(n,k));
    
    if nsamp > maxSamp
        nsamp = maxSamp;
    end
    
    if mod(k,500) == 0
        fprintf(1,'CN: k = %i \n',k);
    end
    
	subsets = f_NcombK(n,k,nsamp);
	nc = size(subsets,1);
	logdk = 0;
	for c = 1:nc
		I = subsets(c,:);
		ld = logdet(COV(I,I));
		if isinf(ld)
			fprintf(2,'ERROR (CN): bad covariance matrix in(k=%d)\n',k);
			err = 2;
			return;
		end
		logdk = logdk + ld;
	end
	cn = cn + (logdk/nc-(k/n)*logd);
end
cn = cn/2;
