clear
close all
clc
    
tic

%% Parameters definition

alpha = 1; %1;
sigma = 1; 
dt = 0.2;  

W  = load('W.dat'); 
N = sqrt(size(W,1));

%% Numerical integration of the dynamical equations

T = 20000; % number of time points

x = zeros(N^2,T);
x(:,1) = randn(N^2,1)*dt;
noise = randn(N^2,T);

for t = 1:T-1
    
    x(:,t+1) = x(:,t) - alpha*x(:,t)*dt + W*x(:,t)*dt + dt*sigma*noise(:,t); 
   
end

toc

%% Calculating principal components

C = cov(x');

[u,s] = eig(C);
s = abs(diag(s));
[s,ind] = sort(s);
u = u(:,ind);
s = flipud(s);
u = fliplr(u);
u = real(u);

%% Plotting PCs

figure
for m=1:16
   subplot(4,4,m)
   I = reshape(u(:,m),N,N);
   I = medfilt2(I,[3 3]);
   imagesc(I)
   axis off
   axis equal; axis tight
end
suptitle('principal components')

toc

%% Playing movie
% 
% cmax = max(max(x));
% cmin = min(min(x));
% 
% figure
% for t = 1:1:T-1
%     
%    B = reshape(x(:,t+1),N,N);
%    B = medfilt2(B,[3 3]);
%    imagesc(B,[cmin cmax]);          % plots on the same scale 
%    %imagesc(B);                     % colorscale changes frame by frame
%    title(['t = ' num2str(t)],'fontsize',16)
%    pause(0.05)
%     
% end
