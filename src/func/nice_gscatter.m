function [h] = nice_gscatter(x,y,g,c,varargin)
h = gscatter(x,y,g,c,varargin{:}) ;

sz = h(1).MarkerSize ;

for i = 1:length(h)
    set(h(i), ... 'MarkerFaceColor','flat',...
        'LineWidth',1.5,'MarkerSize',sz*3,...
        'MarkerFaceColor',c(i,:),...
        'MarkerEdgeColor',[0.5 0.5 0.5]) 
end