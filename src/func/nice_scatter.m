function [ho,r,p] = nice_scatter(x,y,varargin)
h = scatter(x,y,...
    varargin{:},...
    'MarkerFaceAlpha',.6,'MarkerFaceColor','flat',...
    'MarkerEdgeAlpha',.9,'LineWidth',1.5) ;
[r,p] = corr(x(:),y(:),'type','s','rows','complete') ;
disp([ 'spearman rho: ' num2str(r) ', pval: ' num2str(p) ])

if nargout > 0
   ho = h ; 
end