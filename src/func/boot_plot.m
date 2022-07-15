function [] = boot_plot(x,y,bootIter,extend)

if nargin < 3
   bootIter = 1000 ;  
end

if nargin < 4
   extend = 0 ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add ones for intercept
X = [x ones(size(x,1),1)] ;

%% full regression
[fullreg,~,~,~,stats] = regress(y,X) ; 

disp(['R^2 is: ' num2str(stats(1))])
disp(['p-val is: ' num2str(stats(3))])

disp(['slope: ' num2str(fullreg(1))])
disp(['y-intercept: ' num2str(fullreg(2))])

xMin = min(x) ; xMax = max(x) ;

if extend
    xLim = xlim ;
    xMin = xLim(1) ; xMax = xLim(2) ;
end

xVec = xMin:0.05:xMax; % important var for drawing

%% now the bootstrapping

% get the slope and intercept 
boot = bootstrp(bootIter, @regress, y, X);

% get the confidence intervals
bootCI = zeros(size(boot,1),length(xVec)) ;
for idx = 1:bootIter
    bootCI(idx,:) = polyval(boot(idx,:),xVec);  
end
p95 = prctile(bootCI,[2.5,97.5]);    

%% add to the scatter (from above)

hold on;
fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.7,'edgealpha',0);

% and the line
yhat = polyval(fullreg,xVec); 
plot(xVec,yhat,'Color',[0.5 0.5 0.5],'linewidth',2.5);

hold off