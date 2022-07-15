% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

filename = [ DD_PROC '/MamiNum_corrected_proc_2.csv' ] ;
ssheet = readtable(filename,'Delimiter','\t') ;

%% read species

filename = [ DD_RAW '/MamiNum_corrected_sciname.csv' ] ;
scisheet = readtable(filename,'Delimiter',',') ;

lil_sci_names = lower(strrep(strtrim(scisheet.scientific_name),' ','_')) ;

n_animal = length(lil_sci_names) ;

%% read in an example tree to get mammal names

tree_dir = [ DD_RAW '/Completed_5911sp_topoCons_NDexp/' ] ;
first_tree = [ tree_dir '/' 'MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre' ] ;

phy1 = phytreeread(first_tree) ;
big_sci_names = lower(strtrim(get(phy1,'LeafNames'))) ;

%% and match 

[~,findInBig]= ismember(lil_sci_names,big_sci_names) ;

%% now gather all 10,000 matricies 

% outMat = nan(n_animal,n_animal,1e4) ;

% %% distances
% 
% [pd,c] = pdist(phy1,'squareform',true) ; 
% 
% smallerpd = pd(findInBig,findInBig) ;

parpool(10)

%% loop it

o_dir = [ '/mnt/drive2/treeout/' ] ;
mkdir(o_dir) 

parsave = @(fname, tr1,tr2) save(fname, 'tr1','tr2','-v7.3') ;

parfor idx = 1:1e4 
   
    disp(idx)
    
    try
    
        ff = [ tree_dir 'MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree' ...
            num2str(idx-1,'%04.f') '.tre' ] ;
        ff2 = [ o_dir '/tree' num2str(idx-1,'%04.f') '.mat' ] ;
        
        if isfile(ff2)
           disp(['already done: ' num2str(idx)]) 
           continue
        end

        pp = phytreeread(ff) ;
        pd1 = pdist(pp,'squareform',1) ; 
        pd2 = pdist(pp,'squareform',1,'Criteria','levels') ; 

        tr1 = pd1(findInBig,findInBig) ;
        tr2 = pd2(findInBig,findInBig) ; 
%         clear pd ;

%         save(ff2,'tr','-v7.3') 
         parsave(ff2,tr1,tr2)
    catch
       warning([ 'SOMETHING wrong for:' num2str(idx) ]) 
    end
    
end

delete(gcp)

%% now gather into single mats 

n_animal = length(lil_sci_names) ;
mm = logical(triu(ones(n_animal),1)) ;

pd1 = zeros(sum(mm,'all'),1e4) ;
pd2 = zeros(sum(mm,'all'),1e4) ;

for idx = 1:1e4 
   
    disp(idx)

    ff2 = [ o_dir '/tree' num2str(idx-1,'%04.f') '.mat' ] ;
        
    ll = load(ff2) ;
       
    pd1(:,idx) = ll.tr1(mm) ;
    pd2(:,idx) = ll.tr2(mm) ;
        
end

%% viz it

for idx = 1:1000
    subplot(1,2,1)
    imagesc(make_square(pd1(:,idx),1))
    title(num2str(idx))
    subplot(1,2,2)
    imagesc(make_square(pd2(:,idx),1))
        title(num2str(idx))
    waitforbuttonpress
end

% for idx = 1:1000
%     imagesc(make_square(pd2(:,idx),1))
%     title(num2str(idx))
%     waitforbuttonpress
% end

%% read in outlier info

filename = [ DD_INTERM '/' OUTSTR '_divergences_n_outlierdetect.mat' ] ;
divo = load(filename) ;

thr_vals = [ 0 0.05 0.1 0.15 ] ; 

%% threshold based on outliers

for tdx = 1:4

    % the animals that passed outlier detection at this thr
    pass_ani = find(~divo.outliers_thr_nonan(:,tdx)) ; % indexes of the original 225
    % number of animals at this level
    n_animal = length(pass_ani) ;

%     nmask = logical(triu(ones(n_animal),1)) ;
    
    phydist1 = zeros(n_animal,n_animal,1e4) ;
    phydist2 = zeros(n_animal,n_animal,1e4) ;
        
    for idx = 1:1e4
    
        disp(idx)
        
        mm = make_square(pd1(:,idx),1) ;
        ss = mm(pass_ani,pass_ani) ;
        phydist1(:,:,idx) = ss ; 
 
        mm = make_square(pd2(:,idx),1) ;
        ss = mm(pass_ani,pass_ani) ;
        phydist2(:,:,idx) = ss ; 
        
    end

    % and now save it
    filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    save(filename,'phydist1','phydist2','-v7.3')
    
end
       
%% add some names (after the fact)

for tdx = 1:4
   
    disp(tdx)
    
    filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    load(filename) 
    
    % the animals that passed outlier detection at this thr
    pass_ani = find(~divo.outliers_thr_nonan(:,tdx)) ; % indexes of the original 225
    % number of animals at this level
    n_animal = length(pass_ani) ;
    
    scinames = lil_sci_names(pass_ani) ;
    
    save(filename,'phydist1','phydist2','scinames','-v7.3')

end

%% one more tree to load

% readtr = phytreeread([ DD_RAW '/just_sci_names.nwk']) ;
% tmp_tr_names = get(readtr,'leafnames') ;
% % matchup the names
% [aa,findInTrOL]= ismember(lil_sci_names,lower(tmp_tr_names)) ;
% findInTrOL(findInTrOL==0) = [] ;
% 
% trdist = pdist(readtr,'squareform',1) ;
% trdist = trdist(findInTrOL,findInTrOL) ;
% 
% trimmedssheet = ssheet(aa,:) ;
% 
% trimdivmat = divo.animal_div_mat(:,:,2) ;
% trimdivmat = trimdivmat(aa,aa) ;
% 
% scatter(triuvec(trimdivmat),triuvec(trdist))
% corr(triuvec(trimdivmat),triuvec(trdist),'type','s','rows','complete')

% %% do a lil lookin?
% 
% for tdx = 1:4
%     
%     filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
%     pp = load(filename) ;
%     
%     selVec = ~divo.outliers_thr_nonan(:,tdx ) ;
% 
%     divMat = divo.animal_div_mat(:,:,tdx) ;
%     dMat = divMat(selVec,selVec) ;
%     
%     %% all points
% %     mm = logical(triu(ones(sum(selVec)),1)) ;
% %     cc = zeros(1e4,1) ;
% %     for idx = 1:1e4
% %         
% %         disp(idx)
% %         
% %         tmp = pp.phydist1(:,:,idx) ;
% %         cc(idx) = corr(dMat(mm),tmp(mm),'type','s') ;
% %         
% %     end
%     
%     %% block mats
%     
%     filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
%     ll = load(filename,'newsheet') ; 
%     ss = ll.newsheet ;
%     gg = grp2idx(ss.Order) ;
%     
%     [~,dMatBL] = get_block_mat(dMat,gg) ;    
% 
%     % block nusiance regressors
%     snrMat = ss.SNR * ss.SNR' ;
%     [~,snrBL] = get_block_mat(snrMat,gg) ;
%     snrBL(isinf(snrBL)) = NaN ;
%     
%     brVolMat = ss.log10_BrV_ * ss.log10_BrV_' ; 
%     [~,brVolBL] = get_block_mat(brVolMat,gg) ;
%     brVolBL(isinf(brVolBL)) = NaN ;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % run it
%     
%     % mask out singleton comparison too
%     mm2 = logical(triu(ones(size(dMatBL)),0)) ;
%     mm2(isnan(sum(dMatBL)),:) = 0 ;
%     mm2(:,isnan(sum(dMatBL))) = 0 ;
%     
%     cc2 = zeros(1e4,1) ;
%     pc2 = zeros(1e4,1) ;
% 
%     aa = zeros(sum(mm2,'a')*1e4,1) ;
%     bb = zeros(sum(mm2,'a')*1e4,1) ;
%     
%     for idx = 1:1e4
%         
%         disp(idx)
%         
%         tmp = pp.phydist1(:,:,idx) ;
%         [~,tmp2] = get_block_mat(tmp,gg) ;
%         
%         cc2(idx) = corr(dMatBL(mm2),tmp2(mm2),'type','s','rows','complete') ;
%         pc2(idx) = partialcorr(dMatBL(mm2),tmp2(mm2),...
%             [ snrBL(mm2) brVolBL(mm2)],...
%             'type','s','rows','complete') ;
% 
%         % add to vecs for big mat
%         ii_1=((idx-1)*sum(mm2,'a'))+1 ; 
%         ii_2 = ii_1 + sum(mm2,'a') -1 ;  
%         
%         aa(ii_1:ii_2) = dMatBL(mm2) ;
%         bb(ii_1:ii_2) = tmp2(mm2) ;
%         
%     end
% 
% end
