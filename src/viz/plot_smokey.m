function [] = plot_smokey(xvals,yvals,ystd,color1,color2,plotpattern)

if nargin < 4
    color1 = [0.5 0.5 0.5] ;
end

if nargin < 5
    color2 = [0.75 0.75 0.75] ;
end

if nargin < 6
   plotpattern = [ 1 1 ] ; 
end

if ~any(plotpattern)
   error('need to plot something yo') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yvals = yvals(:)' ;
ystd = ystd(:)' ;
xvals = xvals(:)' ; 

if plotpattern(2)
    fill([xvals fliplr(xvals)],[yvals + ystd fliplr(yvals - ystd)],...
        color2,'facealpha',.05,'edgealpha',0);
end

if plotpattern(1)
    hold on
    plot(xvals,yvals,'Color',[ color1 0.5],'linewidth',2.5)
    hold off
end
