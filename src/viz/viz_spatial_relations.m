% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% setup vars

thr_vals = [ 0 0.05 0.1 0.15 ] ; 
odir = [ PROJ_DIR '/reports/figures/commspace/' ] ;
mkdir(odir)

for tDx = 1:4

    %% some viz stuff

    labs = { 'Efficency' '-Searchinfo' 'Diff. efficency' 'Communicability' } ; 

    filename = [ DD_PROC '/' OUTSTR '_distcorr_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    load(filename)
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tDx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    ssheet = ll.newsheet ;

    %% setup

    ff = [ PROJ_DIR '/reports/lab_meeting/animal_colors.mat' ] ;
    load(ff,'anicmap')

    %% loop it

    comm_names = { 'Shortestpaths' 'NegSearchInfo' 'DiffusionEff' 'Communicability' } ;

    for idx = 1:4

        dotsize = 100 ;
        
        f = figure(...
        'units','inches',...
        'position',[0 0 12 6],...
        'paperpositionmode','auto');

        subplot(1,2,1)

    % % %     xx = [] ;
    % % %     yy = [] ;
    % % %     
    % % %     for jdx = 1:nReps
    % % %         
    % % %         xx = [ xx ; ssheet.log10_BrV_] ;
    % % %         yy = [ yy ; squeeze(distcorr_wei_null(:,idx,jdx))] ;
    % % %     end
    % % %     
    % % %     rr = randsample(length(xx),1e4,0) ;
    % % %     
    % % %     nice_scatter(xx(rr),yy(rr),200,repmat([0.9 0.9 0.9],length(xx(rr)),1)) 
    % % %     hold on 
    % % %     nice_scatter(ssheet.log10_BrV_,distcorr_wei(:,idx),200,grp2idx(ssheet.Order))     
    % % %     hold off 

        yy = [] ;

        nReps=size(distcorr_wei_null,3) ;
        for jdx = 1:nReps
            yy = [ yy squeeze(distcorr_wei_null(:,idx,jdx))] ; %#ok<AGROW>
        end    

        maxyy = prctile(yy,2.5,2) ;
        minyy = prctile(yy,97.5,2) ; 

        % just the size
        xx = ssheet.log10_BrV_ ;

        plotx = [] ;
        ploty = [] ;
        for ii = 1:length(minyy)

            tx = [ xx(ii) xx(ii) NaN ] ;
            ty = [ minyy(ii) maxyy(ii) NaN ] ;
            plotx = [ tx, plotx] ; %#ok<AGROW>
            ploty = [ ty, ploty ] ;  %#ok<AGROW>

        end

        plot(plotx,ploty,'Color',[0.9 0.9 0.9],'LineWidth',3) ;

        hold on     
        nice_scatter(ssheet.log10_BrV_,distcorr_wei(:,idx),dotsize,anicmap(grp2idx(ssheet.Order),:))
        hold off 

        xl = xlim ;
        xrange = xl(2) - xl(1) ;

        xlim([ xl(1)-(xrange*0.01) xl(2)+(xrange*0.01)])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subplot(1,2,2)

        xx = distcorr_wei(:,idx) ;
    %     
    %     nice_scatter(ssheet.log10_BrV_,xx,200,grp2idx(ssheet.Order)) 
    %     title(['brain volume vs. ' labs{idx} '~euclidean'])
    %     waitforbuttonpress

        % and plot the corr versus permuted corr
        rr = corr(ssheet.log10_BrV_,xx,'type','s','rows','complete') ;

        nReps = size(distcorr_wei_null,3) ;
        nd = nan(nReps,1) ;
        for jdx = 1:nReps
            nd(jdx) = corr(ssheet.log10_BrV_,...
                squeeze(distcorr_wei_null(:,idx,jdx)),'type','s',...
                'rows','complete') ;
        end
        histogram(nd,'FaceColor',[0.9 0.9 0.9],'LineStyle','none',...
            'Normalization','probability')
        hold on 
        title(['emp. ' comm_names{idx} '~euclidean vs null'])
        
        yl = ylim ;
        xl = xlim ;
        
        line([rr, rr], [ 0 yl(2)*.95], 'LineWidth', 2, 'Color', [1 0.55 0]);
                

        pp = (sum(nd<rr)+1)/(nReps+1) ;
        disp([ 'pval: ' num2str(pp) ])

        yrange = yl(2)-yl(1) ;
        xrange = xl(2) - xl(1) ;
        
        text(rr+(xrange*0.05),(yl(2)*.9),{[ '\rho: ' num2str(round(rr,4)) ]})
        text(rr+(xrange*0.05),(yl(2)*.9)-(yrange*0.05),{[ '{\itp}: ' num2str(round(pp,4)) ]})

        hold off
        

%         waitforbuttonpress

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save it
        
        ff = [ odir '/sp_vs' comm_names{idx} '_thr' num2str(thr_vals(tDx)) '.pdf' ] ;
        print(gcf(),'-dpdf',ff,'-bestfit');
        close all

    end

end % tdx

%% 

% % nice_scatter(ssheet.log10_BrV_,distcorr_sp,200,grp2idx(ssheet.Order)) 
% % title('brain volume vs. shortestpaths~euclidean')
% % 
% % nice_scatter(ssheet.log10_BrV_,distcorr_si,200,grp2idx(ssheet.Order)) 
% % title('brain volume vs. -searchinfo~euclidean')
% % 
% % nice_scatter(ssheet.log10_BrV_,distcorr_mt,200,grp2idx(ssheet.Order)) 
% % title('brain volume vs. diffeff~euclidean')
% % 
% % nice_scatter(ssheet.log10_BrV_,distcorr_co,200,grp2idx(ssheet.Order)) 
% % title('brain volume vs. communicability~euclidean')
% 
% % same results with 
% nice_scatter(ssheet.log10_BrV_,distcorr_bin(:,1),200,grp2idx(ssheet.Order)) 
% nice_scatter(ssheet.log10_BrV_,distcorr_bin(:,2),200,grp2idx(ssheet.Order)) 
% nice_scatter(ssheet.log10_BrV_,distcorr_bin(:,3),200,grp2idx(ssheet.Order)) 
% nice_scatter(ssheet.log10_BrV_,distcorr_bin(:,4),200,grp2idx(ssheet.Order)) 
% 



