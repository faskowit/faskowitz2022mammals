% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% setup vars

thr_vals = [ 0 0.05 0.1 0.15 ] ; 

ff = [ PROJ_DIR '/reports/lab_meeting/animal_colors.mat' ] ;
load(ff)

odir = [ PROJ_DIR '/reports/figures/netdist/' ] ;
mkdir(odir)

for tDx = 1:4

%% some viz stuff
    
    filename = [ DD_PROC '/' OUTSTR '_repanimal_rewrdist_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    load(filename)
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    ssheet = ll.newsheet ;

    %%

    dotsize = 100 ;
    
    f = figure(...
        'units','inches',...
        'position',[0,0,12,6],...
        'paperpositionmode','auto');
    
    ts = tight_subplot(1,3,0.02,[0.1 0.1]) ;
    
%     subplot(1,3,1) 
    axes(ts(1))
    [~,r,p] = nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),geo_dist ),...
        dotsize,anicmap(grp2idx(ssheet.Order),:)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 
    ylabel('Network portrait divergence')
    text(1.75, 0.975,[ '\rho: ' num2str(round(r,4)) ]);
    text(1.75, 0.95,[ '{\itp}: ' num2str(round(p,4)) ]);
    title('Space & topology preserve')

%     subplot(1,3,2) 
    axes(ts(2))
    [~,r,p] = nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),rand_dist ),...
        dotsize,anicmap(grp2idx(ssheet.Order),:)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 
    set(gca,'YTick',[])
    xlabel('Log10 brain volume')
    text(1.75, 0.975,[ '\rho: ' num2str(round(r,4)) ]);
    text(1.75, 0.95,[ '{\itp}: ' num2str(round(p,4)) ]);
    title('Topology preserve')

%     subplot(1,3,3) 
    axes(ts(3))
    [~,r,p] = nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),randmio_dist ),...
        dotsize,anicmap(grp2idx(ssheet.Order),:)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 
    set(gca,'YTick',[])
    text(1.75, 0.975,[ '\rho: ' num2str(round(r,4)) ]);
    text(1.75, 0.95,[ '{\itp}: ' num2str(round(p,4)) ]);
    title('Degree preserve')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'thr' num2str(thr_vals(tDx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff,'-bestfit');
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% for idx = 2:21
% 
%     nice_scatter(...
%         cellfun(@(x_)x_(idx,1),rand_dist ),...
%         cellfun(@(x_)x_(idx,1),geo_dist ),200,grp2idx(ssheet.Order)) 
% %         nice_scatter(...
% %             cellfun(@(x_)x_(idx,1),rand_dist ),...
% %             cellfun(@(x_)x_(idx,1),geo_dist ),200,ones(n_animal,1).*idx) 
%     title([ 'randomize: ' num2str(randlevels(idx)) ])
%     xlim([0 0.5])
%     ylim([0 0.5])
%     hold on 
%     refline(1,0)
% %         hold off 
% 
%     xlabel('random rewire')
%     ylabel('geo rewire')
% 
% %         waitforbuttonpress
% 
% end



%%

aa = grp2idx(ssheet.Order) ;
n_ani_classes = length(unique(grp2idx(ssheet.Order))) ;
% ani_colors = parula(n_ani_classes) ;
ani_names = unique(ssheet.Order) ;

for idx = 1:n_ani_classes

    disp(idx)
    
    tmpcell = geo_dist(aa==sIdx2(idx)) ;
    
    subplot(3,4,sIdx2(idx)) 
    for jdx = 1:size(tmpcell,1) 

        col = anicmap(sIdx2(idx),:) ;
        
        hold on
        % function [] = plot_smokey(xvals,yvals,ystd,color1,color2,plotpattern)
        plot_smokey(randlevels,tmpcell{jdx}(:,1),tmpcell{jdx}(:,2),...
            col,col,[1 1])
        ylim([ 0 0.4])        
    end
            title(ani_names{sIdx2(idx)})
    
end

hold off
