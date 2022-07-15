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

odir = [ PROJ_DIR '/reports/figures/phydist/' ] ;
mkdir(odir)

triunroll = @(x_) x_(logical(triu(ones(size(x_,1)),1))) ;

distlist = {'bin-gen-max' 'bin-gen-mean' 'wei-gen-max' 'wei-gen-mean' ...
    'lap-spec-full' 'lap-spec-5' 'lap-spec-50' 'lap-spec-100' ...
    'lapspec-js' ...
    'adj-spec-full' 'adj-spec-5' 'adj-spec-50' 'adj-spec-100' ...
    'wei-netsimile' 'bin-netsimile' 'net-pd' ...
    } ;

%% loop over threshold

for tdx = 1:4
    
%     filename = [ DD_PROC '/phydistcorr_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
%     if isfile(filename)
%         disp('already done')
%        continue 
%     end
    
    % load the phy distacnes
    filename = [ DD_PROC '/phydist_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    pp = load(filename) ;
    
    
    selVec = ~divo.outliers_thr_nonan(:,tdx ) ;
    divMat = divo.animal_div_mat(:,:,tdx) ;
    netdMat = divMat(selVec,selVec) ;
    
    trilmask = logical(tril(ones(size(netdMat,1)))) ;
%     triumask = logical(tril(ones(size(netdMat,1)),1)) ;
    netdMat(trilmask) = 0 ;
    netdMat = triu(netdMat,1) + triu(netdMat,1)' ;  
    
    % load the other distances
    filename = [ DD_PROC '/' OUTSTR '_altdistances_thr' num2str(THRDENS(tdx)) '_.mat' ] ;
    load(filename) % this loads distMats
  
    % add the netdMat to the other distMats
    distMats(size(distMats,2)+1).mat = netdMat ; 
        
    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    ll = load(filename,'newsheet') ; 
    ss = ll.newsheet ;
    gg = grp2idx(ss.Order) ;
    nBL = max(gg) ;
    
    % also load spatial embedding parameter
    filename = [ DD_PROC '/' OUTSTR '_generative2_thr' num2str(thr_vals(tdx)) '_.mat' ] ;
    gen = load(filename) ;
    
    n_animal = size(netdMat,1) ;
    
    %% all points
        
    utmask = logical(triu(ones(sum(selVec)),1)) ;
    
    %% block mats
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % block nusiance regressors
    % snrMat = ss.SNR * ss.SNR' ;
    snrMat = zeros(n_animal) ;
    for idx = 1:n_animal ; for jdx = 1:n_animal %#ok<ALIGN>
        if idx > jdx ; continue ; end
            snrMat(idx,jdx) = sqrt(ss.SNR(idx)*ss.SNR(jdx)) ;
    end ; end
    snrMat = triu(snrMat,1) + triu(snrMat,1)' ;
   
    [~,snrBL] = get_block_mat(snrMat,gg) ;
    snrBL(isinf(snrBL)) = NaN ;
    
%     brVolMat = ss.log10_BrV_ * ss.log10_BrV_' ; 
    brVolMat = zeros(n_animal) ;
    for idx = 1:n_animal ; for jdx = 1:n_animal %#ok<ALIGN>
        if idx > jdx ; continue ; end
            brVolMat(idx,jdx) = sqrt(ss.log10_BrV_(idx)*ss.log10_BrV_(jdx)) ;
    end ; end
    brVolMat = triu(brVolMat,1) + triu(brVolMat,1)' ;

    [~,brVolBL] = get_block_mat(brVolMat,gg) ;
    brVolBL(isinf(brVolBL)) = NaN ;
    
%     genPar = gen.estPar * gen.estPar' ;
    genPar = zeros(n_animal) ;
    for idx = 1:n_animal ; for jdx = 1:n_animal %#ok<ALIGN>
        if idx > jdx ; continue ; end
            genPar(idx,jdx) = sqrt(gen.estPar(idx)*gen.estPar(jdx)) ;
    end ; end
    genPar = triu(genPar,1) + triu(genPar,1)' ;
    
    [~,parBL] = get_block_mat(genPar,gg) ;
    parBL(isinf(parBL)) = NaN ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,dMatBL] = get_block_mat(netdMat,gg) ;    
    nDist = size(distMats,2) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run it
    
%     nreps = 1e3 ;
    nreps = 1e4 ; 
    
    % mask out singleton comparison too
    blmask = logical(triu(ones(nBL),0)) ;
    blmask(isnan(sum(dMatBL)),:) = 0 ;
    blmask(:,isnan(sum(dMatBL))) = 0 ;
    
    FULL_rhovals = zeros(nreps,nDist) ;
    FULL_prhovals_1 = zeros(nreps,nDist) ;
    FULL_prhovals_2 = zeros(nreps,nDist) ;
    
    BL_rhovals = zeros(nreps,nDist) ;
    BL_prhovals_1 = zeros(nreps,nDist) ;
    BL_prhovals_2 = zeros(nreps,nDist) ;

    aa = zeros(sum(blmask,'a')*nreps,nDist) ;
    bb = zeros(sum(blmask,'a')*nreps,nDist) ;
    
    for idx = 1:nreps
        
        disp(idx)
        
        tmp = pp.phydist1(:,:,idx) ;
        [~,tmp2] = get_block_mat(tmp,gg) ;
        
        for jdx = 1:nDist
        
            % full 
            
            ddd = distMats(jdx).mat ;
            
            % full
            FULL_rhovals(idx,jdx) = corr(ddd(utmask),tmp(utmask),'type','s','rows','complete') ;
            FULL_prhovals_1(idx,jdx) = partialcorr(ddd(utmask),tmp(utmask),...
                [ snrMat(utmask) brVolMat(utmask)],...
                'type','s','rows','complete') ;
            FULL_prhovals_2(idx,jdx) = partialcorr(ddd(utmask),tmp(utmask),...
                [ snrMat(utmask) genPar(utmask)],...
                'type','s','rows','complete') ; 
            
            % block
            
            [~,dBL] = get_block_mat(distMats(jdx).mat,gg) ; 
            
            BL_rhovals(idx,jdx) = corr(dBL(blmask),tmp2(blmask),'type','s','rows','complete') ;
            BL_prhovals_1(idx,jdx) = partialcorr(dBL(blmask),tmp2(blmask),...
                [ snrBL(blmask) brVolBL(blmask)],...
                'type','s','rows','complete') ;
            BL_prhovals_2(idx,jdx) = partialcorr(dBL(blmask),tmp2(blmask),...
                [ snrBL(blmask) parBL(blmask)],...
                'type','s','rows','complete') ;

            % add to vecs for big mat
            ii_1=((idx-1)*sum(blmask,'a'))+1 ; 
            ii_2 = ii_1 + sum(blmask,'a') -1 ;  

            aa(ii_1:ii_2,jdx) = dBL(blmask) ;
            bb(ii_1:ii_2,jdx) = tmp2(blmask) ;
        
        end
        
    end

    %% and now save it
    filename = [ DD_PROC '/phydistcorr_thr' num2str(thr_vals(tdx)) '.mat'  ] ;
    save(filename,'distlist','distMats','BL_*','FULL_*','aa','bb','-v7.3')
    
end

%% look at it

% for idx = 1:nDist
%    histogram(BL_prhovals_2(:,idx))
%    title(distlist{idx})
%    waitforbuttonpress
%     
% end
% 
% 
% for idx = 1:nDist
%    histogram(BL_prhovals_2(:,idx))
%    hold on
%     
% end



