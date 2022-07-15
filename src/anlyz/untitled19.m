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

%% load centoid data from previosly

% filename = [ DD_INTERM '/' OUTSTR '_divergences_n_outlierdetect.mat' ] ;
% netportdivs = load(filename,'animal_div_mat*');

% filename = [ DD_INTERM '/' OUTSTR '_just_cent_ani_dat.mat' ] ;
% centdata = load(filename,'leastdist*','verydist*',...
%               'outliers_thr','outliers_thr_nonan') ;

%% get the data

centmat = struct() ;
for idx = 1:4

    filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(THRDENS(idx)) '_.mat' ] ;
    ll = load(filename) ; 
    centmat(idx).mat = ll.data ;
    centmat(idx).ssheet = ll.newsheet ;
          
end
clear ll ;              
          
%% get the distances

centdist = struct() ;
for idx = 1:4

    filename = [ DD_INTERM '/cents_repani_stack_thr' num2str(THRDENS(idx)) '_.mat' ] ;
    ll = load(filename) ; 
    centdist(idx).dist = ll.data ;
          
end
clear ll ;

%% make distances 

% test it
amat = centmat(1).mat(:,:,1) ;
bmat = centmat(1).mat(:,:,2) ;

amatd = squareform(pdist(centdist(1).dist(:,:,1))) ;
bmatd = squareform(pdist(centdist(1).dist(:,:,2))) ;

% test some distances

[dmax,dmean,K] = genmodel_dist_bin(amat,bmat,amatd,bmatd) ;
[d2max,d2mean,K2] = genmodel_dist_wei(amat,bmat,amatd,bmatd) ;

lapdist = f_lap_dist(amat,bmat) ;

lapdist2 = f_lap_dist(amat,bmat,5) ;
lapdist3 = f_lap_dist(amat,bmat,10) ;
lapdist4 = f_lap_dist(amat,bmat,20) ;

jsdist = f_JS_discr(amat,bmat) ;

%%

tdx = 2 ; 
n_animal = size(centmat(tdx).mat,3) ;

cMat = centmat(tdx).mat ;
dMat = centdist(tdx).dist ;

% distMats = struct() ;
% for ddx = 1:9
%    distMats(ddx).mat = nan(n_animal) ;  
% end
% distMats = cell(9,1) ;
% for ddx = 1:9
%    distMats{ddx} = zeros(n_animal) ;  
% end

mat1 = nan(n_animal) ;
mat2 = nan(n_animal) ;
mat3 = nan(n_animal) ;
mat4 = nan(n_animal) ;
mat5 = nan(n_animal) ;
mat6 = nan(n_animal) ;
mat7 = nan(n_animal) ;
mat8 = nan(n_animal) ;
mat9 = nan(n_animal) ;

parfor idx = 1:n_animal
    
    amat = cMat(:,:,idx) ;
    amatd = squareform(pdist(dMat(:,:,idx))) ;

    % clean up mat
    gc = get_components(amat) ;
    gcmode = mode(gc) ;
    if any(gc~=gcmode(1))
        % make a lil smaller 
        amat = amat(gc==gcmode(1),gc==gcmode(1)) ;
        amatd = amatd(gc==gcmode(1),gc==gcmode(1)) ; 
    end     
    
    resVec = nan(n_animal,9) ;
    
    for jdx = 1:n_animal
        if idx >= jdx ; continue ; end
        
        disp([ num2str(idx) ' - ' num2str(jdx) ]) 
                
        bmat = cMat(:,:,jdx) ;    
        bmatd = squareform(pdist(dMat(:,:,jdx))) ;

        % clean up mat
        gc = get_components(bmat) ;
        gcmode = mode(gc) ;
        if any(gc~=gcmode(1))
            % make a lil smaller 
            bmat = bmat(gc==gcmode(1),gc==gcmode(1)) ;
            bmatd = bmatd(gc==gcmode(1),gc==gcmode(1)) ; 
        end    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % binary genmodel dist
        [d1max,d1mean] = genmodel_dist_bin(amat,bmat,amatd,bmatd) ;
        mat1(idx,jdx) = d1max ;
        mat2(idx,jdx)  = d1mean ;

        % weighted genmodel dist
        [d2max,d2mean] = genmodel_dist_wei(amat,bmat,amatd,bmatd) ;
        mat3(idx,jdx)  = d2max ;
        mat4(idx,jdx)  = d2mean ;
        
        % lap dist
        mat5(idx,jdx)  = f_lap_dist(amat,bmat) ;
        
        % lap dist n components
        mat6(idx,jdx)  = f_lap_dist(amat,bmat,5) ;
        mat7(idx,jdx)  = f_lap_dist(amat,bmat,10) ;
        mat8(idx,jdx)  = f_lap_dist(amat,bmat,20) ;

        % lap js div
        mat9(idx,jdx)  = f_JS_discr(amat,bmat) ;
        
    end
end

%%

mm = logical(tril(ones(n_animal),0)) ;

distMats = struct() ;
distMats(1).mat = mat1 ; 
distMats(2).mat = mat2 ; 
distMats(3).mat = mat3 ; 
distMats(4).mat = mat4 ; 
distMats(5).mat = mat5 ; 
distMats(6).mat = mat6 ; 
distMats(7).mat = mat7 ; 
distMats(8).mat = mat8 ; 
distMats(9).mat = mat9 ; 

for ddx = 1:9
   tmp = distMats(ddx).mat ;
   tmp(mm) = 0 ;
   distMats(ddx).mat = tmp + tmp' ; 
end

% mat1(mm) = 0 ;
% mat1 = mat1 + mat1' ; 

for idx = 1:9

    imsc_grid_comm(distMats(idx).mat,grp2idx(centmat(tdx).ssheet.Order))
    waitforbuttonpress

end
