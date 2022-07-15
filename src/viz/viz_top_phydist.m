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

odir = [ PROJ_DIR '/reports/figures/top_phy/' ] ;
mkdir(odir)

triunroll = @(x_) x_(logical(triu(ones(size(x_,1)),1))) ;

distlist = {'bin-gen-max' 'bin-gen-mean' 'wei-gen-max' 'wei-gen-mean' ...
    'lap-spec-full' 'lap-spec-5' 'lap-spec-50' 'lap-spec-100' ...
    'lapspec-js' ...
    'adj-spec-full' 'adj-spec-5' 'adj-spec-50' 'adj-spec-100' ...
    'wei-netsimile' 'bin-netsimile' 'net-pd' ...
    } ;

rhos_names = { 'rho0' 'prho1' 'prho2' } ;
rhos_longnames = { 'Rank corr. (\rho)' ...
    'Partial rank corr. (snr + brain vol.)' ...
    'Partial rank corr. (snr + gen. param.)' } ;

topN = 0.01 * 10e3 ;

%% loop over threshold

for tdx = 1:4
    
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    load(filename,'newsheet') ; 
    ssheet = newsheet ; 
    
    % load the phy distacnes
    filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    pp = load(filename) ;
    [u,ia,ic] = unique(pp.scinames) ;
    
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
    
%     % load the other distances
%     filename = [ DD_PROC '/' OUTSTR '_altdistances_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
%     load(filename)
%   
%     % add the netdMat to the other distMats
%     distMats(size(distMats,2)+1).mat = netdMat ; 

    n_animal = size(netdMat,1) ;

    %% the phydistances 
    
    filename = [ DD_PROC '/phydistcorr_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    load(filename)

    %% 
    
    f_resstruct = struct() ;
    bl_resstruct = struct() ;
    
    f_resstruct(1).vals = FULL_rhovals ;
    f_resstruct(2).vals = FULL_prhovals_1 ;
    f_resstruct(3).vals = FULL_prhovals_2 ;
   
    bl_resstruct(1).vals = BL_rhovals ;
    bl_resstruct(2).vals = BL_prhovals_1 ;
    bl_resstruct(3).vals = BL_prhovals_2 ;
    
    %% block

    for fdx = 1:3
    
        analzrho = bl_resstruct(fdx).vals(:,16) ;
        % let's just look at netpd
        [sortedRho,sss] = sort(bl_resstruct(fdx).vals(:,16),'descend') ;
        
        topMinMax = [ min(sortedRho(1:topN)) max(sortedRho(1:topN)) ] ;
        botMinMax = [ min(sortedRho((end-topN+1):end)) ...
                      max(sortedRho((end-topN+1):end)) ] ;
        
        % get the consensus from the top trees
        toptrees = pp.phydist1(:,:,sss(1:topN)) ;
        bottrees = pp.phydist1(:,:,sss((end-topN+1):end)) ;

%         toptrees = mean(pp.phydist1(:,:,sss(1:topN)),3) ;
%         bottrees = mean(pp.phydist1(:,:,sss((end-topN+1):end)),3) ;
        
%         toptrees = triu(toptrees,1) + triu(toptrees,1)' ; 
%         bottrees = triu(bottrees,1) + triu(bottrees,1)' ; 
        
        top_modz = zeros(n_animal,topN) ;
        bot_modz = zeros(n_animal,topN) ;
%         top_v = [] ;
%         bot_v = [] ;
%         n_order = max(spord.vec) ;
%         top_bm = zeros((n_order*(n_order-1)/2 + n_order),topN) ;
%         bot_bm = zeros((n_order*(n_order-1)/2 + n_order),topN) ;

        for idx = 1:topN
            
            tt = toptrees(:,:,idx) ;
            tt = tt + tt' ;
            bb = bottrees(:,:,idx) ;
            bb = bb + bb' ;
            
            top_modz(:,idx) = module_degree_zscore(tt,spord.vec) ; 
            bot_modz(:,idx) = module_degree_zscore(bb,spord.vec) ;
            
%             top_v = [triunroll(tt) ; top_v ];
%             bot_v = [triunroll(bb) ; bot_v ];
%             
%             [~,top_bl] = get_block_mat(tt,spord.vec) ;
%             triunroll(top_bl)
            
        end
        
%         scatter(triunroll(mean(toptrees,3)),triunroll(mean(bottrees,3)))
        
%         histscatter(prct_discretize(triunroll(toptrees),100),...
%                 prct_discretize(triunroll(bottrees),100),100)
        f = figure(...
            'units','inches',...
            'position',[0,0,10,8],...
            'paperpositionmode','auto');

        ts = tight_subplot(2,3,0.04,0.02,0.1) 
        
        axes(ts(1))
        
        top_d = mean(toptrees,3) ;
        top_d = top_d + top_d' ;
        [~,reorm ] = reorder_matrix(top_d,'line',500) ;
%         reorm = 1:190 ;

        [h,~,nt] = imsc_grid_comm(top_d(reorm,reorm),spord.vec(reorm),[],[1 0.55 0],0,spord.cat) ;
        axis square
        caxis([0 200])
        colormap(brewermap(100,'YlGnBu'))
        set(gca,'xtick',[])
%         set(gca,'ytick',[])
        cb = colorbar()
        cb.Visible = 'off'

        title({'Top corr. trees' , [ '(' num2str(round(topMinMax(1),2)) ...
               '-' num2str(round(topMinMax(2),2)) ')' ] })

        
        axes(ts(2))
        
        bot_d = mean(bottrees,3) ;
        bot_d = bot_d + bot_d' ;
%         [~,reorm ] = reorder_matrix(top_d,'circ',2000) ;
        
        [h,~,nt] = imsc_grid_comm(bot_d(reorm,reorm),spord.vec(reorm),[],[1 0.55 0],0) ;
        axis square 
        caxis([0 200]) 
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap(brewermap(100,'YlGnBu'))
        cb = colorbar()

        title({'Bottom corr. trees' , [ '(' num2str(round(botMinMax(1),2)) ...
               '-' num2str(round(botMinMax(2),2)) ')' ] })    
           
           
        axes(ts(3))
        
%         h = scatter(triunroll(top_d),triunroll(bot_d),'filled',...
%             'MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
%         xlim([0 350]) ; ylim([0 350])
%         refline(1,0)
%         cb = colorbar()
%         cb.Visible = 'off' 
        h = histscatter(triunroll(top_d),triunroll(bot_d),30)
        cm = brewermap(100,'blues') ; 
        cm2 = [ 1 1 1 ; cm(50:100,:) ] ;
        colormap(ts(3),cm2)
        axis square 
        caxis([0 25])
        cb = colorbar()
        
        axis square
        grid minor
        
        [r,p] = corr(triunroll(top_d),triunroll(bot_d),'type','s') ;
      
        text(0.94,0.03,[ '\rho: ' num2str(round(r,3)) ...
            ', \it{p}' '-value: ' num2str(round(p,3)) ], ...
            'Units','normalized',...
            'Rotation',90)
        
%         histscatter(prct_discretize(top_v,100),prct_discretize(bot_v,100),20)
%         histscatter(top_v,bot_v,100)
%
%         tt = mean(toptrees,3) ;
%         tt = tt + tt' ; 
%         oo = doDiffusionMap(corr(tt),n_animal,2) ; 
%         
%         tt = mean(bottrees,3) ;
%         tt = tt + tt' ; 
%         ii = doDiffusionMap(corr(tt),n_animal,2) ;    
          

%     [~,bm] = get_block_mat(dd,grp2idx(ssheet.Order)) ;

        axes(ts(4))

         [~,bm] = get_block_mat(top_d,spord.vec) ;

        nnn = isnan(sum(bm,2)) ;
        bm = bm(~nnn,~nnn) ; 

    %     bm(isnan(bm)) = 0 ; 

        % subplot(1,2,2)
        imsc_grid_comm(bm,1:(max(spord.vec)-sum(nnn)),3,[1 0.55 0])
        cb = colorbar() ; cb.Visible = 'off' ;
        colormap(ts(4),brewermap(100,'YlGnBu'))
        caxis([0 200])
        axis square
        set(gca,...
            'xtick',1:(max(spord.vec)-sum(nnn)),...
            'xticklabel',spord.cat(~nnn));
        set(gca,...
            'ytick',1:(max(spord.vec)-sum(nnn))) 
        
        axes(ts(5))

         [~,bm] = get_block_mat(bot_d,spord.vec) ;

        nnn = isnan(sum(bm,2)) ;
        bm = bm(~nnn,~nnn) ; 

    %     bm(isnan(bm)) = 0 ; 

        % subplot(1,2,2)
        imsc_grid_comm(bm,1:(max(spord.vec)-sum(nnn)),3,[1 0.55 0])
        colorbar
        colormap(ts(5),brewermap(100,'YlGnBu'))
        caxis([0 200])
        axis square
        set(gca,...
            'xtick',1:(max(spord.vec)-sum(nnn)),...
            'xticklabel',spord.cat(~nnn));
        set(gca,...
            'ytick',1:(max(spord.vec)-sum(nnn))) 
        
        
        axes(ts(6))
        
        h = histscatter(top_modz(:),bot_modz(:),50)
        colormap(ts(6),flipud(brewermap(100,'BuPu')))
        axis square 
        caxis([0 400])
        cb = colorbar()

        h.ShowEmptyBins = 'off'
        
        xlabel('Top corr. trees')
        ylabel('Bottom corr. trees')
        
        grid minor
        
        title('Mod. deg. z-score')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save the figure
        set(0, 'DefaultFigureRenderer', 'painters');
        ff = [ odir 'topbot_block_' rhos_names{fdx} '_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
        print(gcf(),'-dpdf',ff,'-bestfit');
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end % fdx


    
end