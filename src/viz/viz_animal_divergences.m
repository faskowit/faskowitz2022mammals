% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% viz between-subjects distances
% 
% filename = [ DD_INTERM '/species_treesim.mat' ] ;
% spsim = load(filename) ; 

%% read in outlier info

filename = [ DD_INTERM '/run1_divergences_n_outlierdetect.mat' ] ;
divo = load(filename) ;

%% setup vars

thr_vals = [ 0 0.05 0.1 0.15 ] ; 
odir = [ PROJ_DIR '/reports/figures/divergences/' ] ;
mkdir(odir)

%%  divergence mats and boxplots 

for tDx = 1:4

    thrIdx = tDx ;     
    
    %% somehow make a new animal order

    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(thrIdx)) '_.mat' ] ;
    load(filename,'newsheet') ; 
    ssheet = newsheet ; 
    
    filename = [ DD_INTERM '/trim_tree_thr' num2str(thr_vals(thrIdx)) '.mat' ]  ;
    if ~isfile(filename)
    
        % load the phy distacnes
        filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(thrIdx)) '.mat'  ] ;
        pp = load(filename,'scinames') ;
        scinames = pp.scinames ;     

        % read in consensus tree
        cons_tr = phytreeread([ DD_RAW '/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre' ]) ;
        cons_tr_names = get(cons_tr,'LeafNames') ;
        cons_tr_names = lower(regexprep(strtrim(cons_tr_names),'[A-Z]{2,}',''));
        cons_tr_names = strrep(cons_tr_names,'__','') ;

        %% manipulate info from tree

        [uniq_scinames,ia,ic] = unique(scinames) ;

        [~,findInBig]= ismember(uniq_scinames,cons_tr_names) ;

        bigpdist = pdist(cons_tr,'Squareform',true) ;
        smallerpdist = bigpdist(findInBig,findInBig) ;

        % reduce the tree   
        trim_tr = seqlinkage(smallerpdist,'weighted',uniq_scinames) ;
        trimtr_leafnames = get(trim_tr,'LeafNames') ;

        [~,trimtrInSci]= ismember(trimtr_leafnames,scinames) ;
        trimtrOrd = ssheet.Order(trimtrInSci) ;

        order_order = unique(trimtrOrd,'stable') ;

        % save the new tree
        filename = [ DD_INTERM '/trim_tree_thr' num2str(thr_vals(thrIdx)) '.nwk' ]  ;
        phytreewrite(filename,trim_tr)
    
        filename = [ DD_INTERM '/trim_tree_thr' num2str(thr_vals(thrIdx)) '.mat' ]  ;
        save(filename,'trim_tr','trimtrOrd')
        
        %% make the species order

        spc_order_vec = nan(size(ssheet.Order)) ; 
    %     spc_order_catvec = cell(size(ssheet.Order)) ;
        for idx = 1:length(order_order) 
            inds = strcmp(ssheet.Order,order_order{idx}) ;
            spc_order_vec(inds) = idx ;
    %         spc_order_catvec(inds) = repmat(cellstr(order_order{idx}),sum(inds),1) ;
        end

        spord = struct() ;
        spord.vec = spc_order_vec ;
        spord.cat = order_order ;
        spord.speciesvec =  ssheet.Order ;

        filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(thrIdx)) '.mat' ] ;
        save(filename,'spord') 
        
    end % if trim_tree exists
  
    filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(thrIdx)) '.mat' ] ;
    load(filename)
       
    filename = [ DD_INTERM '/trim_tree_thr' num2str(thr_vals(thrIdx)) '.mat' ]  ;
    ttt = load(filename) ;
    
    %% viz it
    
    f = figure(...
        'units','inches',...
        'position',[0,0,12,6],...
        'paperpositionmode','auto');

%     [trg,trgg] = grp2idx(trimtrOrd) ;
    mmm = pdist(ttt.trim_tr,'Squareform',true') ;
    [h,~,nt] = imsc_grid_comm(mmm,grp2idx(ttt.trimtrOrd),[],[1 0.55 0],0,spord.cat) ;

    colormap()
    %set(gca,'xtick',[])
%     set(gca,'YTickLabelRotation',30)
    set(gca,'xtick',nt,'XTickLabel',spord.cat) 

    colormap(brewermap(100,'YlGnBu'))
    cb = colorbar() ;
    cb.Label.String = 'Patristic distance' ;
%     cb.Ticks = [] ;
    axis square

%     axes(ts(2)) 
% 
% %     mmm2 = ranktrnsf_und(mmm) ;
%     mmm2 = mmm .* 1;
%     mmm2(1:(size(mmm,1)+1):end) = 0 ;
%     [~,sortidx ] = sort(ssheet.BrVol(trimtrInSci)) ;
%     imagesc(mmm2(sortidx,sortidx))
% %     caxis([0 200])
%     colorbar
%     axis square
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save the figure
    ff = [ odir 'phydist_thr' num2str(thr_vals(thrIdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %% make the species order
% 
%     spc_order_vec = nan(size(ssheet.Order)) ; 
% %     spc_order_catvec = cell(size(ssheet.Order)) ;
%     for idx = 1:length(order_order) 
%         inds = strcmp(ssheet.Order,order_order{idx}) ;
%         spc_order_vec(inds) = idx ;
% %         spc_order_catvec(inds) = repmat(cellstr(order_order{idx}),sum(inds),1) ;
%     end
% 
%     spord = struct() ;
%     spord.vec = spc_order_vec ;
%     spord.cat = order_order ;
%     spord.speciesvec =  ssheet.Order ;
%     
%     filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(thrIdx)) '.mat' ] ;
%     save(filename,'spord') 
    
    %%
    
%     filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(thrIdx)) '_.mat' ] ;
%     ll = load(filename) ; 
%     ssheet = ll.newsheet ;
%     dat = ll.data ;

    selVec = ~divo.outliers_thr_nonan(:,thrIdx ) ;

    % ss = spsim.speciesSim(selVec,selVec) ;
    % imsc_grid_comm(ss,grp2idx(ssheet.Order))

    divMat = divo.animal_div_mat(:,:,thrIdx) ;

    dd = divMat(selVec,selVec) ;
    dd(~~tril(ones(size(dd,1)),-1)) = 0;
    dd = dd + dd' ;
    dd(1:size(dd,1)+1:end) = 0 ;
        
    % subplot(1,2,1)
    [d,di] = reorderMAT(-dd,1e3,'line') ;
%     [g,gg] = grp2idx(ssheet.Order) ;
    [h,~,nt] = imsc_grid_comm(dd(di,di),spord.vec(di),3,[1 0.55 0],0,spord.cat) ;
    colorbar
    colormap(brewermap(100,'YlGnBu'))
    caxis([ 0 0.5 ])
    axis square
    set(gca,'xtick',nt,'XTickLabel',spord.cat) 
    % set(gca,'ytick',[])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'btwn_ani_mat_thr' num2str(thr_vals(thrIdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot block
    
%     [~,bm] = get_block_mat(dd,grp2idx(ssheet.Order)) ;
     [~,bm] = get_block_mat(dd,spord.vec) ;

    nnn = isnan(sum(bm,2)) ;
    bm = bm(~nnn,~nnn) ; 
    
%     bm(isnan(bm)) = 0 ; 

    % subplot(1,2,2)
    imsc_grid_comm(bm,1:(max(spord.vec)-sum(nnn)),3,[1 0.55 0])
    colorbar
    colormap(brewermap(100,'YlGnBu'))
    caxis([ 0 0.5 ])
    axis square
    set(gca,...
        'xtick',1:(max(spord.vec)-sum(nnn)),...
        'xticklabel',spord.cat(~nnn));

%     set(gca,'xtick',[]) ; 
%     set(gca,'ytick',[])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'btwn_ani_block_thr' num2str(thr_vals(thrIdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    onD = bm(~~eye(size(bm,1))) ;
    offD = triuvec(bm,1) ;

    figure
    cc = brewermap(4,'Paired') ;
    fcn_boxpts([onD ; offD], [ ones(length(onD),1) ; ones(length(offD),1)*2 ], ...
        cc(3:4,:),0,{ 'Within' 'Between' })
    ylim([0 0.5])

    % and do the between test
    p = ranksum(onD,offD,'tail','left') ;

    text(1.025,.35,...
        ['\it{U}-test \it{p}: ' num2str(round(p,3)) ], ...
        'Rotation',90,...
        'Units','normalized')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'btwn_boxplot_thr' num2str(thr_vals(thrIdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % tDx

%% sort orders into size colors

filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(3)) '_.mat' ] ;
ll = load(filename) ; 
ssheet = ll.newsheet ;
filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(3)) '.mat' ] ;
load(filename) ;

% dat = ll.data ;
% % make density really quick
% dens = nan(size(dat,3),1) ;
% for idx = 1:size(dat,3)
%     dens(idx) = density_und(dat(:,:,idx)) ;
% end
% [g,gg] = grp2idx(ssheet.Order) ;

sz_by_group = parcellate_lab(ssheet.log10_BrV_,spord.vec) ;
[ss,sIdx] = sort(sz_by_group,'ascend') ; 
[ss2,sortAniGroups] = sort(sIdx) ; 

ac = flipud(brewermap(max(spord.vec),'spectral')) ;

anicmap = ac(sortAniGroups,:) ; 

fcn_boxpts(ssheet.log10_BrV_,spord.vec,anicmap,0,spord.cat) 
ylabel('Log_{10} Brain Volume')

ff = [ PROJ_DIR '/reports/figures/order_sizes.pdf' ] ;
print(gcf(),'-dpdf',ff);
close all

ff = [ PROJ_DIR '/reports/figures/animal_colors.mat' ] ;
save(ff,'anicmap','sortAniGroups')

%% do the gen parameter plots

for tdx = 1:4
    
    filename = [ DD_PROC '/' OUTSTR '_generative2_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    gen = load(filename) ;
    
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    
    filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(tdx)) '.mat' ] ;
    load(filename) ;
    
%     nice_scatter(ll.newsheet.log10_BrV_,gen.estPar,200,spord.vec) 
%     colormap(anicmap)
%     

    
    f = figure(...
        'units','inches',...
        'position',[0,0,8,8],...
        'paperpositionmode','auto');
    
    hold on
    for idx = 1:max(spord.vec)
%         s(idx) = scatter(0,0,0.1,anicmap(idx,:),'filled') ;
        nice_scatter(ll.newsheet.log10_BrV_(spord.vec==idx),...
            gen.estPar(spord.vec==idx),200,anicmap(idx,:)) 

    end
      
    [r,p] = corr(ll.newsheet.log10_BrV_,gen.estPar,'type','s') ;
   
    
    text(0.7,0.92,[ '\rho: ' num2str(round(r,3)) ...
        ', \it{p}' '-value: ' num2str(round(p,3)) ], ...
        'Units','normalized')
    
    lgd = legend(spord.cat,'NumColumns',2,...
        'Location','southwest') ;
    
    ylabel('Spatial parameter')
    xlabel('Log_{10} Brain Volume (mm^3)')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'param_scatter_' num2str(thr_vals(tdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
end

%% now do the extended ttest figure

for tdx = 1:4
    
%     % load the phy distacnes
%     filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
%     pp = load(filename) ;
  
    filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(tdx)) '.mat' ] ;
    load(filename,'spord') 

    selVec = ~divo.outliers_thr_nonan(:,tdx ) ;
    divMat = divo.animal_div_mat(:,:,tdx) ;
    netdMat = divMat(selVec,selVec) ;
    
    trilmask = logical(tril(ones(size(netdMat,1)))) ;
%     triumask = logical(tril(ones(size(netdMat,1)),1)) ;
    netdMat(trilmask) = 0 ;
    netdMat = triu(netdMat,1) + triu(netdMat,1)' ;  
    
    % load the other distances
    filename = [ DD_PROC '/' OUTSTR '_altdistances_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    load(filename)
  
    % add the netdMat to the other distMats
    distMats(size(distMats,2)+1).mat = netdMat ; 
    
    %%
    
    triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ;
    
    ndist = size(distMats,2) ;
    distcolors = brewermap(ndist,'greens') ;
    
    withinvec = [] ;
    btwnvec = [] ;
    
    for idx = 1:ndist

        rr = ranktrnsf_und(distMats(idx).mat) ;
        
        [~,bm] = get_block_mat(rr,spord.vec) ;

        nnn = isnan(sum(bm,2)) ;
        bm = bm(~nnn,~nnn) ; 
        
        withinvec = [ withinvec(:) ; diag(bm) ] ;
        btwnvec = [ btwnvec(:) ; triunroll(bm) ] ; 
        
    end
    
    cc = brewermap(4,'Paired') ;
    fcn_boxpts([withinvec ; btwnvec], [ ones(length(withinvec),1) ; ones(length(btwnvec),1)*2 ], ...
    cc(3:4,:),0,{ 'Within' 'Between' })

    ylim([0 18e3 ])

    ylabel('Mean block rank')

    % and do the between test
    p = ranksum(withinvec,btwnvec,'tail','left') ;

    text(1.025,.35,...
        ['\it{U}-test \it{p}: ' num2str(round(p,3)) ], ...
        'Rotation',90,...
        'Units','normalized')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the figure
    ff = [ odir 'alldist_btwn_boxplot_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
    print(gcf(),'-dpdf',ff);
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    close all
    
end


