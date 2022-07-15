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

%% read in outlier info

filename = [ DD_INTERM '/' OUTSTR '_divergences_n_outlierdetect.mat' ] ;
divo = load(filename) ;

%% compute the masks for each thr

thr_vals = [ 0 0.05 0.1 0.15 ] ; 

filename = [ DD_INTERM '/con_mat_gn_kperm_stack.mat' ] ;
ll = load(filename) ;

%%

thr_masks = struct() ;
na = size(ll.data,3) ;

for tdx = 2:4
   
    pass_ani = find(~divo.outliers_thr_nonan(:,tdx)) ;
    n_animal = length(pass_ani) ;
    % get the centroid for each animal we are keeping
    leastd = divo.leastdist(:,tdx) ;
    
    thr_masks(tdx).mat = zeros(NNODES,NNODES,n_animal) ;
    
    for idx = 1:n_animal
        
        disp(idx)
        
        curr_a_ind = pass_ani(idx) ;
        curr_least = leastd(curr_a_ind) ; 
        
        thr_masks(tdx).mat(:,:,idx) = ...
            fcn_mst_plus(ll.data(:,:,curr_a_ind,curr_least),thr_vals(tdx),'inv') ;
        
    end
    
end

%% noooooowww save some new data
% the data that needs to be thresholded based on imposed density... 

thr_val_str = { '0' '0.05' '0.1' '0.15' } ;

load_dats = { 'con_mat_gn_kperm_stack.mat' ...
    'con_mat_kperm_stack.mat' 'leng_mat_kperm_stack.mat' } ;

for idx = 1:length(load_dats)
    
    % load it
    filename = [ DD_INTERM '/' load_dats{idx} ] ;
    ll = load(filename) ;
    ddd = ll.data ; 
    
    for tdx = 1:length(thr_vals)
        
        disp([ num2str(idx) ' ' num2str(tdx)]) 
        
        out_name = [ strrep(load_dats{idx},'kperm_stack.mat','') ...
        'repani_stack_thr' thr_val_str{tdx} '_.mat' ] ;
        
        % the animals that passed outlier detection at this thr
        pass_ani = find(~divo.outliers_thr_nonan(:,tdx)) ; % indexes of the original 225
        % number of animals at this level
        n_animal = length(pass_ani) ;
        % get the centroid for each animal 
        leastd = divo.leastdist(:,tdx) ; % length 225
        
        data = nan(size(ddd,1),size(ddd,2),n_animal) ;
        
        if tdx > 1 % then we should apply threshold
           
            
            for kdx = 1:n_animal
                disp(kdx)
                curr_a_ind = pass_ani(kdx) ;
                curr_least = leastd(curr_a_ind) ;  
                
                data(:,:,kdx) = single(thr_masks(tdx).mat(:,:,kdx)) .* ddd(:,:,curr_a_ind,curr_least) ; 
            end
            
        else
            for kdx = 1:n_animal
                disp(kdx)
                curr_a_ind = pass_ani(kdx) ;
                curr_least = leastd(curr_a_ind) ;  
                
                data(:,:,kdx) = ddd(:,:,curr_a_ind,curr_least) ; 
            end
            
        end
        
        newsheet = ssheet(pass_ani,:) ;
        
        save(out_name,'data','newsheet','-v7.3') 
        
    end
    
end

%% and the other data

thr_val_str = { '0' '0.05' '0.1' '0.15' } ;

load_dats = { 'cents_kperm_stack.mat' 'node_adj_kperm_stack.mat' } ;

for idx = 1:length(load_dats)
    
    % load it
    filename = [ DD_INTERM '/' load_dats{idx} ] ;
    ll = load(filename) ;
    ddd = ll.data ; 
    
    for tdx = 1:length(thr_vals)
        
        disp([ num2str(idx) ' ' num2str(tdx)]) 
        
        out_name = [ strrep(load_dats{idx},'kperm_stack.mat','') ...
        'repani_stack_thr' thr_val_str{tdx} '_.mat' ] ;
        
        % get the outliers
        outlr = divo.outliers_thr(:,tdx) ;
        keep = ~outlr ;
    
        % get the centroid for each animal we are keeping
        ld = divo.leastdist(keep,tdx) ;
        
        newddd = ddd(:,:,keep,ld) ;
        
        save(out_name,newddd,'-v7.3') 
        
    end
    
end




