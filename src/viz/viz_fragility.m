% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% setup

ff = [ PROJ_DIR '/reports/figures/animal_colors.mat' ] ;
load(ff,'anicmap')

%% setup vars

thr_vals = [ 0 0.05 0.1 0.15 ] ; 
odir = [ PROJ_DIR '/reports/figures/fragility/' ] ;
mkdir(odir)

for tDx = 1:4

    dotsize = 100 ;

    f = figure(...
    'units','inches',...
    'position',[0 0 8 6],...
    'paperpositionmode','auto');
    
    %% some viz stuff

    filename = [ DD_PROC '/' OUTSTR '_fragilityMyDef2_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    load(filename)
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    ssheet = ll.newsheet ;

    %% some viz stuff

    fragstack = zeros(30,length(hub_frag)) ;
    for idx = 1:length(hub_frag)
       fragstack(:,idx) = hub_frag{idx}(~isnan(hub_frag{idx})) ;
    end

    [~,r,p] = nice_scatter(ssheet.log10_BrV_,mean(fragstack),dotsize,anicmap(grp2idx(ssheet.Order),:)) 
    xl = xlim() ;
    xrange = xl(2) - xl(1) ;
    yl = ylim() ;
    yrange = yl(2) - yl(1) ;
    
    text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.1),[ '\rho: ' num2str(round(r,4)) ]);
    text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.05),[ '{\itp}: ' num2str(round(p,4)) ]);
    
    % save the figure
    ff = [ odir '/fragility_thr' num2str(thr_vals(tDx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all

    % nice_scatter(ssheet.log10_BrV_,std(hfstack),200,grp2idx(ssheet.Order)) 
%     [a,b] = corr(ssheet.log10_BrV_,mean(fragstack)'); 

    f = figure(...
    'units','inches',...
    'position',[0 0 8 6],...
    'paperpositionmode','auto');

    % entropy 
    [~,r,p] = nice_scatter(ssheet.log10_BrV_,hub_ent_trap,dotsize,anicmap(grp2idx(ssheet.Order),:)) 
    xl = xlim() ;
    xrange = xl(2) - xl(1) ;
    yl = ylim() ;
    yrange = yl(2) - yl(1) ;
    
    text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.1),[ '\rho: ' num2str(round(r,4)) ]);
    text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.05),[ '{\itp}: ' num2str(round(p,4)) ]);
    
    % save the figure
    ff = [ odir '/frag_ent_thr' num2str(thr_vals(tDx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all

end

%% stack all fragilities

allfrag = [] ;
allbrvol = [] ;

for tDx = 1:4
    
    filename = [ DD_PROC '/' OUTSTR '_fragilityMyDef2_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    load(filename)
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    ssheet = ll.newsheet ;

    %% some viz stuff

    fragstack = zeros(30,length(hub_frag)) ;
    for idx = 1:length(hub_frag)
       fragstack(:,idx) = hub_frag{idx}(~isnan(hub_frag{idx})) ;
    end    
    
    mf = mean(fragstack) ;
    allfrag = [ allfrag ; mf(:) ] ;
    allbrvol = [ allbrvol ; ssheet.log10_BrV_ ] ;

end

h = histscatter(allbrvol,allfrag,[30 30])
axis square
h.LineStyle = 'none' ;
h.ShowEmptyBins = 'on' ;
cb = colorbar
colormap(brewermap(100,'BuPu'))
cb.Label.String = 'Density of points' ;
xlabel('Log10 Brain Volume')
ylabel('Mean fragility')

[r,p] = corr(allbrvol,allfrag,'type','s')

xl = xlim() ;
xrange = xl(2) - xl(1) ;
yl = ylim() ;
yrange = yl(2) - yl(1) ;

text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.1),[ '\rho: ' num2str(round(r,4)) ]);
text(xl(1)+(xrange*0.03), yl(1)+(yrange*0.05),[ '{\itp}: ' num2str(round(p,5)) ]);

% save the figure
ff = [ odir '/all_histscatter.pdf' ] ;
print(gcf(),'-dpdf',ff);
close all

%%

% filename = [ DD_PROC '/' OUTSTR '_repanimal_hubresil_thr' num2str(thr_vals(thrIdx)) '_.mat' ] ;
% load(filename)
% 
% resilstack = zeros(200,length(hub_resil)) ;
% for idx = 1:length(hub_resil)
%    resilstack(:,idx) = hub_resil{idx}(~isnan(hub_resil{idx})) ;
% end
% 
% %%
% 
% nice_scatter(ssheet.log10_BrV_,median(resilstack),200,anicmap(grp2idx(ssheet.Order),:)) 
% 
% % nice_scatter(mean(fragstack),median(resilstack),200,ssheet.log10_BrV_) 
% 




