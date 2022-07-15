% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% setup vars

thr_vals = [ 0 0.05 0.1 0.15 ] ; 

filename = [ DD_INTERM '/' OUTSTR '_divergences_n_outlierdetect.mat' ] ;
divo = load(filename) ;

ff = [ PROJ_DIR '/reports/lab_meeting/animal_colors.mat' ] ;
load(ff)

odir = [ PROJ_DIR '/reports/figures/outliers/' ] ;
mkdir(odir)

triunroll = @(x_) x_(logical(triu(ones(size(x_,1)),1))) ;

topN = 0.01 * 10e3 ;

%% loop over threshold

for tdx = 1:4
    
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    load(filename,'newsheet') ; 
    ssheet = newsheet ; 
    
%     % load the phy distacnes
%     filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
%     pp = load(filename) ;
%     [u,ia,ic] = unique(pp.scinames) ;
    
    filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(tdx)) '.mat' ] ;
    load(filename,'spord') 

    ff = [ PROJ_DIR '/reports/figures/animal_colors.mat' ] ;
    load(ff,'anicmap','sortAniGroups')
 
    selVec = ~divo.outliers_thr_nonan(:,tdx ) ;
    divMat = divo.animal_div_mat(:,:,tdx) ;
    netdMat = divMat(selVec,selVec) ;
    
    trilmask = logical(tril(ones(size(netdMat,1)))) ;
%     triumask = logical(tril(ones(size(netdMat,1)),1)) ;
    netdMat(trilmask) = 0 ;
    netdMat = triu(netdMat,1) + triu(netdMat,1)' ;  

    n_animal = size(netdMat,1) ;
    
    %%
  
    f = figure(...
        'units','inches',...
        'position',[0,0,10,6],...
        'paperpositionmode','auto');
    
    divMat = triu(divMat,1) + triu(divMat,1)' ;
    
    ss = ~selVec+1 ; 
%     [rr] = reorder_mod(divMat,ss)
    
    [~,rr] = reorderMAT(divMat(rr,rr),100,'line') ;
    
    strrr = cell(2,1) ;
    strrr{1} = [ 'dataset (' num2str(sum(selVec)) ')' ] ; 
    strrr{2} = [ 'outliers (' num2str(sum(~selVec)) ')' ] ; 
    [h,iii] = imsc_grid_comm(divMat(rr,rr),ss(rr),3,[],[],strrr) ;
    set(h,'AlphaData',~isnan(divMat(iii,iii)))
    
%     colormap(flipud(brewermap(100,'Blues')))
    colormap(brewermap(100,'Blues'))
    colorbar
        
    set(gca,'Color',[0.5 0.5 0.5])
    
    caxis([0 1])
    
    axis square
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save the figure
    ff = [ odir 'outl_mat_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    
    aa = nonzeros(divMat(selVec,selVec)) ;
    bb = divMat(selVec,~selVec) ;
    
%     cc = brewermap(2,'Dark2') ;
    
%     fcn_boxpts([aa(:) ; bb(:)], [ ones(length(aa(:)),1) ; ones(length(bb(:)),1)*2 ], ...
%         cc(1:2,:),0,{ 'Within dataset' 'Between dataset & outliers' })
%     ylim([0 0.5])
    
    histogram(aa(:)) 
    hold on
    histogram(bb(:))
   
    % and do the between test
    p = ranksum(aa(:),bb(:)) ;

    text(1.025,.35,...
        ['\it{U}-test \it{p}: ' num2str(round(p,3)) ], ...
        'Rotation',90,...
        'Units','normalized')
    
    
    legend({ 'Within dataset' 'Between dataset and outliers' })
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save the figure
    ff = [ odir 'outl_hist_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end