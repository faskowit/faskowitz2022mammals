% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% read in proc spreadsheet

filename = [ DD_PROC '/MamiNum_corrected_proc.txt' ] ;
ssheet = readtable(filename,'Delimiter','\t') ;

animal_names = regexprep(ssheet.file_prefix,'_FF3_allfibers_kperm','') ;
n_animal = length(animal_names) ;

VIZIT = 0 ;

%% load data

filename = [ DD_INTERM '/con_mat_gn_kperm_stack.mat' ] ;
ll = load(filename) ;
allCon = ll.data ; clear ll ;
disp('loaded conn')

filename = [ DD_INTERM '/leng_mat_kperm_stack.mat' ] ;
ll = load(filename) ;
allLen = ll.data ; clear ll ;
disp('loaded lenn')

filename = [ DD_INTERM '/cents_kperm_stack.mat' ] ;
ll = load(filename) ;
allCents = ll.data ; clear ll ;
disp('loaded cent')

% load communbities
filename = [ DD_INTERM '/' OUTSTR '_consensus_comms.mat' ] ;
ll = load(filename,'o_str') ;
commStr = ll.o_str ; clear ll ;

disp('loaded comms')


%% load up the exemplar mats for each animal

filename = [ DD_INTERM '/' OUTSTR '_divergences.mat' ] ;
ll = load(filename) ;
[~,centInd] = min(ll.sort_inds,[],2) ;
clear ll


%% other params

triumask = make_triumask(NNODES) ;
nbins =  100 ;

compare_bins = [ .05 .10 .20 ] ;

out_filename = [ DD_INTERM '/' OUTSTR '_agrVdist.mat' ] ;


if ~isfile(out_filename) 

agr_bin_wei = zeros(nbins,n_animal) ;
agr_max = zeros(n_animal,1) ;
arg_dist_corr = cell(n_animal,1) ;

for idx = 1:n_animal

    tmp_agr = zeros(nbins,NRESAMP) ;
    tmp_agrdistcorr = zeros(NRESAMP,1) ;
    
    aCon = squeeze(allCon(:,:,idx,:)) ;
    aCent = squeeze(allCents(:,:,idx,:)) ;
    aAmat = commStr(idx).Aall ;
    
    parfor jdx = 1:NRESAMP
    
        if mod(jdx,10)==1 ; disp([ num2str(idx) ' - ' num2str(jdx)]) ; end

%             c = aCon(:,:,jdx) ; 
            cents = aCent(:,:,jdx) ;
            d = squareform(pdist(cents)) ;
            amat = aAmat(:,:,jdx) ;
                          
            % get distance bins
%             [~,ee] = histcounts(d(triumask),nbins) ;
%             ee = linspace(min(d(triumask),[],'all'),max(d(triumask),[],'all'),nbins+1) ;
            ee = [-inf quantile(d(triumask),nbins-1)  inf ] ;
            ddisc = discretize(d,ee) ;
            ddisc(1:NNODES+1:end) = NaN ;
            
            % wei in euc bins
            tmp_agr(:,jdx) = parcellate_lab(amat(triumask),ddisc(triumask),(1:nbins)','sum') ;

            tmp_agrdistcorr(jdx) = corr(amat(triumask),d(triumask),'type','s') ;
            
            short_thr = compare_bins
            
%             bin_thr = lower(compare_bins(1)*nbins) ;
%             upper_mas = ddisc>=(nbins-bin_thr) ;
%             lower_mas = ddisc<=bin_thr ;
            
    end % reps
    
    agr_bin_wei(:,idx) = mean(tmp_agr,2) ;
    [~,agr_max(idx)]= max(mean(tmp_agr,2)) ;
    arg_dist_corr{idx} = tmp_agrdistcorr ;

end

else
   
    load(out_filename)
end

%%

nice_scatter(ssheet.log10_BrV_,iqr(agr_bin_wei,1),...
    300,grp2idx(ssheet.Order))
xlabel('log10 brainvol')
ylabel('meas')

nice_scatter(ssheet.log10_BrV_,std(agr_bin_wei,1),...
    300,grp2idx(ssheet.Order))
xlabel('log10 brainvol')
ylabel('meas')

% % scatter(iqr(agr_bin_wei,1),ssheet.log10_BrV_,300,grp2idx(ssheet.Order),'filled') 
% scatter(ssheet.log10_BrV_,agr_max,...
%     300,grp2idx(ssheet.Order),...
%     'MarkerFaceAlpha',.6,'MarkerFaceColor','flat',...
%     'MarkerEdgeAlpha',.9,'LineWidth',1.5)

% scatter(ssheet.log10_BrV_,cellfun(@median,arg_dist_corr),300,grp2idx(ssheet.Order),'filled') 
scatter(ssheet.log10_BrV_,arrayfun(@(x) median(arg_dist_corr{x}) ,1:n_animal),...
    300,grp2idx(ssheet.Order),...
    'MarkerFaceAlpha',.6,'MarkerFaceColor','flat',...
    'MarkerEdgeAlpha',.9,'LineWidth',1.5)

scatter(ssheet.log10_BrV_,arrayfun(@(x) arg_dist_corr{x}(centInd(x)) ,1:n_animal),...
    300,grp2idx(ssheet.Order),...
    'MarkerFaceAlpha',.6,'MarkerFaceColor','flat',...
    'MarkerEdgeAlpha',.9,'LineWidth',1.5)

%% save it

filename = [ DD_INTERM '/' OUTSTR '_agrVdist.mat' ] ;
save(filename,'arg_dist_corr','agr_bin_wei','agr_max','-v7.3')

