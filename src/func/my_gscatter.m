function [h] = my_gscatter(x,y,g,c,varargin)

x = x(:) ;
y = y(:) ;
g = g(:) ;

grps = unique(g) ;

for i = 1:length(grps)
    cg = grps(i) ;
    h(i) = scatter(x(g==cg),y(g==cg),...
            varargin{:},...
            'MarkerFaceColor', c(i,:) ,...
            'MarkerEdgeColor', c(i,:) * 0.9, ...
            'MarkerFaceAlpha',.6,...
            'MarkerEdgeAlpha',.9,...
            'LineWidth',1.5)
    hold on
end
hold off
