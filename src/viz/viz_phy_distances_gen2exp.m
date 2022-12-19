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

odir = [ PROJ_DIR '/reports/figures/phydist_gen2exp/' ] ;
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
    'Partial rank corr. (snr)' ...
    'Partial rank corr. (snr + exp. gen. param.)' } ;

%% loop over threshold

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

    %% now make a picture of all the distance matrices
    ndistmeas = size(distMats,2) ;
    btwndistDist = nan(ndistmeas) ;
    
    for idx = 1:ndistmeas
        for jdx = 1:ndistmeas
            if idx >= jdx ; continue ; end
            
            a = triunroll(distMats(idx).mat) ;
            b = triunroll(distMats(jdx).mat) ;
            
            btwndistDist(idx,jdx) = corr(a,b,'type','Spearman','rows','complete');      
        end
    end
    btwndistDist(1:(ndistmeas+1):end) = 1 ;
 
%     f = figure(...
%         'units','inches',...
%         'position',[0,0,12,6],...
%         'paperpositionmode','auto');
%     % https://www.mathworks.com/matlabcentral/answers/308892-how-to-add-grid-to-the-image
%     h = imagesc(ones(ndistmeas)) ;
%     set(h,'AlphaData',zeros(ndistmeas))
%     [rows, columns] = size(ones(ndistmeas));
%     stepSize = 0.5 ;
%     hold on;
%     lineSpacing = 20; % Whatever you want.
%     for row = stepSize : stepSize : rows
%         line([0, columns+stepSize], [row, row], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
%     end
%     for col = stepSize : stepSize : columns
%         line([col, col], [0, rows+stepSize], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
%     end
%     h = imsc_grid_comm(btwndistDist,1:ndistmeas,0.1,[1 1 1],1,distlist) ;
%     set(h,'AlphaData',~isnan(btwndistDist))
%     colormap(brewermap(100,'BuPu'))
%     cb = colorbar() ;
%     cb.Label.String = "Rank correlation (\rho)" ;
%     caxis([0 1])
%     axis square
%    	set(gca,'xtick',[])
%     hold off
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% save the figure
%     ff = [ odir 'altdist_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
%     print(gcf(),'-dpdf',ff,'-bestfit');
%     close all
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %% the phydistances 
    
    filename = [ DD_PROC '/phydistcorr_extras_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    load(filename)

    %%
    
    % distcolors = brewermap(ndistmeas,'Paired') ;
%     distcolors = brewermap(ndistmeas,'RdYlGn') ;
    % distcolors = brewermap(nDist,'Spectral') ;
    % distcolors = turbo(nDist) ;

    % make categorical colormap
    distcolors = zeros(16,3) ;
    pairedbase = brewermap(12,'Paired') ;
    distcolors(1:4,:) = interp_cmap(pairedbase(1,:),...
                                 pairedbase(2,:), 4) ;
    distcolors(5:9,:) = interp_cmap(pairedbase(3,:),...
                                 pairedbase(4,:), 5) ;                             
    distcolors(10:13,:) = interp_cmap(pairedbase(7,:),...
                                 pairedbase(8,:), 4) ;     
    distcolors(14:15,:) = interp_cmap(pairedbase(5,:),...
                                 pairedbase(6,:), 2) ;     
    distcolors(16,:) = pairedbase(10,:) ;      
    %% 
    
    f_resstruct = struct() ;
    bl_resstruct = struct() ;
    
    f_resstruct(1).vals = FULL_rhovals ;
    f_resstruct(2).vals = FULL_prhovals_1 ;
    f_resstruct(3).vals = FULL_prhovals_2 ;
   
    bl_resstruct(1).vals = BL_rhovals ;
    bl_resstruct(2).vals = BL_prhovals_1 ;
    bl_resstruct(3).vals = BL_prhovals_2 ;
    
    %% full

    for fdx = 1:3
    
        analzrho = f_resstruct(fdx).vals ;

        % calculate the range
        rmin = min(analzrho,[],'all') ;
        rmax = max(analzrho,[],'all') ;
        binsz = (rmax-rmin)./50 ;

        f = figure(...
            'units','inches',...
            'position',[0,0,16,6],...
            'paperpositionmode','auto');
        hold on
        for idx = 1:ndistmeas

            vv = analzrho(:,idx) ;
            pval = prctile(vv,1) ;
            ppass = pval>0 ; 

            h = histogram(vv,'BinWidth',binsz,...
                'Normalization','probability','LineStyle','none') ;
            if ~ppass
               h.FaceColor = [ 0.8 0.8 0.8 ] ;
               h.FaceAlpha = 0.3 ;
            else
               h.FaceColor = distcolors(idx,:) ;
               h.FaceAlpha = 0.5  ;
            end

        end
        legend(distlist,'Location','eastoutside')
        hold off

%         xlim([-0.1 0.3 ])
%         ylim([0 0.2])
        
        ylabel('Normalized count')
        xlabel(rhos_longnames{fdx})
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save the figure
        ff = [ odir 'disthist_full_' rhos_names{fdx} '_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
        print(gcf(),'-dpdf',ff,'-bestfit');
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end % fdx
    
    %% blocks

    for fdx = 1:3
    
        analzrho = bl_resstruct(fdx).vals ;

        rmin = min(analzrho,[],'all') ;
        rmax = max(analzrho,[],'all') ;
        binsz = (rmax-rmin)./50 ;

        f = figure(...
            'units','inches',...
            'position',[0,0,16,6],...
            'paperpositionmode','auto');
        hold on
        for idx = 1:ndistmeas

            vv = analzrho(:,idx) ;
            pval = prctile(vv,1) ;
            ppass = pval>0 ; 

            h = histogram(vv,'BinWidth',binsz,...
                'Normalization','probability','LineStyle','none') ;
            if ~ppass
               h.FaceColor = [ 0.8 0.8 0.8 ] ;
               h.FaceAlpha = 0.3 ;
            else
               h.FaceColor = distcolors(idx,:) ;
               h.FaceAlpha = 0.5  ;
            end

        end
        legend(distlist,'Location','eastoutside')
        hold off

%         xlim([-0.1 0.3 ])
%         ylim([0 0.2])
        
        ylabel('Normalized count')
        xlabel(rhos_longnames{fdx})
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save the figure
        ff = [ odir 'disthist_blocks_' rhos_names{fdx} '_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
        print(gcf(),'-dpdf',ff,'-bestfit');
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    %% look at the histscatter
    
%     f = figure(...
%         'units','inches',...
%         'position',[0,0,10,10],...
%         'paperpositionmode','auto');
% 
%     ts = tight_subplot(4,4,0.05) ;
%     
%     for idx = 1:ndistmeas
%    
%         axes(ts(idx))
%         h =histscatter(prct_discretize(aa(:,idx),100),...
%             prct_discretize(bb(:,idx),100),[20 20]) ;
%         axis square
%         h.ShowEmptyBins = 'on' ; 
%         h.LineStyle = '-' ;
%         h.LineWidth = 0.5 ;
%         colormap([brewermap(100,'YlGn')])
% %         colormap([brewermap(100,'GnBu')])
%         hl = refline(1,0) ;
%         hl.Color = [ 0.3 0.3 0.3 ] ;
%         title(distlist{idx})
%         
%         caxis([ 0 1e4])
%         
%         [r,p] = corr(aa(:,idx),bb(:,idx),'type','p') ;
%         [rho,prho] = corr(aa(:,idx),bb(:,idx),'type','s') ; 
% 
%         if p<1e-6
%             sigp = '*' ;
%         else
%             sigp = '' ;
%         end
%         
%         if prho<1e-6
%             sigph = '*' ;
%         else
%             sigph = '' ;
%         end
%             
%         text(105,5,...
%             ['\it{r}: ' num2str(round(r,3)) sigp ...
%             ', \rho: ' num2str(round(rho,3)) sigph ],...
%             'Rotation',90)
%         
% %         annotation('textbox',[0.9 0.1 0.1 0.1], ...
% %             'string',[ distlist{idx} ' ' num2str(r) ' ' num2str(p) ])
%         
%         if idx<=12
%             set(gca,'xtick',[])
%         end
%         if mod(idx,4)~=1
%             set(gca,'ytick',[])
%         end
%       
%     if idx==16
%         f = figure ;
%         h = imagesc((0:1e4)') ;
%         ylim([0 1e4])
%         colormap(brewermap(100,'YlGn'))
%         set(gca,'xtick',[])
%         set(gca,'ytick',0:1000:1e4)
% %             set(gca,'YTickLabel',num2str((0:1000:1e4)','%1.0e'))
%         axis xy  
%         cb = colorbar ;
%         cb.Label.String = "Bin count" ;
% %             axis off    
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% save the figure
%         ff = [ odir 'histscat_colorbar_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
%         print(gcf(),'-dpdf',ff,'-bestfit');
%         close(f)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           
%     end
%     
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% save the figure
%     ff = [ odir 'histscat_full_thr' num2str(thr_vals(tdx)) '.pdf' ] ;
%     print(gcf(),'-dpdf',ff,'-bestfit');
%     close all
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
end

%% make a plot comparing power vs exp 

for tdx = 1:4

    %%
    
 filename = [ DD_PROC '/' OUTSTR '_generative2_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
gen1 = load(filename) ;

 filename = [ DD_PROC '/' OUTSTR '_generative2exp_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
gen2 = load(filename) ;

filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
 
% filename = [ DD_PROC '/species_order_thr' num2str(thr_vals(tdx)) '.mat' ] ;
% load(filename) ;

ff = [ PROJ_DIR '/reports/figures/animal_colors.mat' ] ;
load(ff,'anicmap','sortAniGroups')

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
    nice_scatter(gen2.estPar(spord.vec==idx),...
        gen1.estPar(spord.vec==idx),200,anicmap(idx,:)) 

end
      
[r,p] = corr(gen1.estPar,gen2.estPar,'type','s') ;


text(0.7,0.92,[ '\rho: ' num2str(round(r,3)) ...
    ', \it{p}' '-value: ' num2str(round(p,3)) ], ...
    'Units','normalized')

lgd = legend(spord.cat,'NumColumns',2,...
    'Location','northwest') ;

xlabel('Exp. spatial parameter')
ylabel('Power spatial parameter')
   
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the figure
ff = [ odir 'scatter_genpar_' num2str(thr_vals(tdx)) '.pdf' ] ;
print(gcf(),'-dpdf',ff);
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


