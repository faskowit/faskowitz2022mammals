function [h,hdl] = fcn_boxpts(x,lab,cmap,samenum,names)
% clear all
% close all
% clc

% load boxdata.mat
% x = mu;

if ~exist('samenum','var')
    samenum = 0;
end
h = hist(lab,1:max(lab));

sz = size(x);
if sz(2) > 1
    lab = repmat(1:sz(2),sz(1),1);
    x = x(:);
    lab = lab(:);
end

if isempty(cmap)
    cmap = jet(max(lab));
end

iqr = [25,50,75];
wid = 0.3;
cmapdark = cmap - 0.25;
cmapdark(cmapdark < 0) = 0;
pointSZ = 100 ;

% f = fcn_rickplot([2,2,4,2]);
% ax = axes;
% hold(ax,'on');
hold on

u = unique(lab);

for j = 1:length(u)
    
    if ~exist('names','var')
        names{j} = sprintf('%i',j);
    end
    
    idx = lab == u(j);
    vals = x(idx);
    prct = prctile(vals,iqr);
    
    whisk1 = min(vals);
    whisk2 = max(vals);
    
    if samenum
        r = randperm(length(vals),min(nonzeros(h)));
        vals = vals(r);
    end
    
    xvals = j + linspace(-wid*0.5,wid*0.5,length(vals));
    s = randperm(length(xvals));
    xvals = xvals(s);
    
    hdl.pointshandle{j} = scatter(xvals,vals,pointSZ,'filled');
    set(hdl.pointshandle{j},...
        'markerfacecolor',cmap(j,:),'markerfacealpha',0.5,...
        'markeredgecolor',cmap(j,:),'markeredgealpha',0.8);
    
    xvals = j + [-wid,+wid,+wid,-wid,-wid];
    yvals = [prct(1),prct(1),prct(3),prct(3),prct(1)];
    hdl.boxhandle{j} = plot(xvals,yvals,'color','k','linewidth',0.75);
    hdl.medianhandle{j} = plot(j + [-wid,wid],prct(2)*ones(1,2),'color','k','linewidth',0.75);
    
    
    xvals = [j*ones(1,2), nan, j*ones(1,2)];
    yvals = [prct(1),whisk1,nan,prct(3),whisk2];
    hdl.whiskerhandle{j} = plot(xvals,yvals,'color','k','linewidth',0.75);
    
end
set(gca,...
    'xlim',[-2*wid + 1,length(u) + 2*wid],...
    'ylim',[min(x) - 0.05*range(x),max(x) + 0.05*range(x)]);
set(gca,'YTickLabelMode','auto')

set(gca,...
    'xtick',1:length(u),'xticklabel',names);
hold off