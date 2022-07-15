% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

VIZIT = 0  %#ok<NOPTS>

%% big loop

for tDx = 1:length(THRDENS)
%%

filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
ll = load(filename) ; 
c_aCon = ll.data ;
ssheet = ll.newsheet ;
clear ll ;

filename = [ DD_INTERM '/leng_mat_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
ll = load(filename) ; 
c_aLeng = ll.data ;
clear ll ;

filename = [ DD_INTERM '/cents_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
ll = load(filename) ; 
c_aCent = ll.data ;
clear ll ;

%% setup

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 
% randlevels = [ 0 ((5:5:100) ./ 100) ] ;
n_animal = size(ssheet,1) ;

m = logical(triu(ones(NNODES),1)) ;

%% init outputs

nReps = 2500 ; 

distcorr_wei = nan(n_animal,4) ;
distcorr_bin = nan(n_animal,4) ;

distcorr_wei_null = nan(n_animal,4,nReps) ;
% distcorr_bin_null = nan(n_animal,4,nReps) ;

%% big loop

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_distcorr_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
if ~isfile(filename)

%%
parfor idx = 1:n_animal

    %%
    
    disp(idx) 
    
    % data 
    c = c_aCon(:,:,idx) ;
    cents = c_aCent(:,:,idx) ;
    eud = squareform(pdist(cents)) ;
    d = c_aLeng(:,:,idx) ;
    
    conmask = c > 0 ;
       
    distcorr_wei(idx,:) = get_comm_corr(c,eud) ;
    distcorr_bin(idx,:) = get_comm_corr(double(c>0),eud) ;

    tic
    for ndx = 1:nReps

        % disp(ndx)

        [~,rewr] = geombinsurr_partial(c,d,1,20,'quantiles') ;

        distcorr_wei_null(idx,:,ndx) = get_comm_corr(rewr,eud) ;
%         distcorr_bin_null(idx,:,ndx) = get_comm_corr(double(rewr>0),eud) ;
    
    end
    toc


    
end

else % if data is already generated, load it
   load(filename) 
end

%% saveit

filename = [ DD_PROC '/' OUTSTR '_distcorr_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'distcorr_*','-v7.3')


%% end the big threshold loop
end

%% some viz stuff

if VIZIT

labs = { 'efficeny' '-searchinfo' 'diffeff' 'communicability' } ; %#ok<UNRCH>
    
thr=0.10;
filename = [ DD_PROC '/' OUTSTR '_distcorr_thr' num2str(thr) '_.mat' ] ;
load(filename)
filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
ssheet = ll.newsheet ;

%% loop it

for idx = 1:4
  
    xx = distcorr_wei(:,idx) ;
    
    nice_scatter(ssheet.log10_BrV_,xx,200,grp2idx(ssheet.Order)) 
    title(['brain volume vs. ' labs{idx} '~euclidean'])
    waitforbuttonpress
    
    % and plot the corr versus permuted corr
    rr = corr(ssheet.log10_BrV_,xx,'type','s','rows','complete') ;
    
    nReps = size(distcorr_wei_null,3) ;
    nd = nan(nReps,1) ;
    for jdx = 1:nReps
        nd(jdx) = corr(ssheet.log10_BrV_,...
            squeeze(distcorr_wei_null(:,idx,jdx)),'type','s',...
            'rows','complete') ;
    end
    histogram(nd)
    hold on 
    title(['emp. ' labs{idx} '~euclidean vs null'])
    line([rr, rr], ylim, 'LineWidth', 2, 'Color', 'r');
    hold off
    
    pp = (sum(nd<rr)+1)/(nReps+1) ;
    disp([ 'pval: ' num2str(pp) ])
    
    waitforbuttonpress
    
end

%% 

nice_scatter(ssheet.log10_BrV_,distcorr_sp,200,grp2idx(ssheet.Order)) 
title('brain volume vs. shortestpaths~euclidean')

nice_scatter(ssheet.log10_BrV_,distcorr_si,200,grp2idx(ssheet.Order)) 
title('brain volume vs. -searchinfo~euclidean')

nice_scatter(ssheet.log10_BrV_,distcorr_mt,200,grp2idx(ssheet.Order)) 
title('brain volume vs. diffeff~euclidean')

nice_scatter(ssheet.log10_BrV_,distcorr_co,200,grp2idx(ssheet.Order)) 
title('brain volume vs. communicability~euclidean')


nice_scatter(ssheet.log10_BrV_,distcorr_sp_bin,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,distcorr_si_bin,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,distcorr_mt_bin,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,distcorr_co_bin,200,grp2idx(ssheet.Order)) 

end


