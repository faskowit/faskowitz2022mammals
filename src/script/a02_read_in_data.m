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
data_dir = [ DD_RAW '/xkperm/' ] ;

%% 

n_animal = size(ssheet,1) ;
all_data = cell(n_animal,1) ;

parfor idx = 1:n_animal
    
    curr = ssheet.file_prefix{idx} ;

    animal_str = struct() ;
    animal_str.name = curr ;
    for jdx = 1:NRESAMP
       
        disp(['animal: ' num2str(idx) ' resamp: ' num2str(jdx) ])
        curr_f = [ data_dir '/' curr num2str(jdx) '.mat' ] ;
        animal_str(jdx).dat = readin_mammal_dat(curr_f) ;  
        
    end
    
    all_data{idx} = animal_str ;
    
end

filename = [ DD_INTERM '/all_data.mat' ] ;
save(filename,'all_data','-v7.3')

%% now write the data into mat files

sqr_entries = { 'dist_mat' 'leng_mat' 'node_adj' 'con_mat' 'con_mat_gn' } ;

for sdx = 1:length(sqr_entries)

    data = zeros(NNODES,NNODES,n_animal,NRESAMP) ;

    % populate data
    for idx = 1:n_animal

        disp(idx)
        
        animal_str = all_data{idx} ;

        for jdx = 1:NRESAMP

            data(:,:,idx,jdx) = animal_str(jdx).dat.(sqr_entries{sdx}) ;

        end
    end
    
    % save data
    filename = [ DD_INTERM '/' sqr_entries{sdx} '_kperm_stack.mat' ] 
    save(filename,'data','-v7.3')
    ls(filename)
   
    clear data
end

%% cents

data = zeros(NNODES,3,n_animal,NRESAMP) ;
for idx = 1:n_animal

    disp(idx)

    animal_str = all_data{idx} ;

    for jdx = 1:NRESAMP

        data(:,:,idx,jdx) = animal_str(jdx).dat.cents ;

    end
end
filename = [ DD_INTERM '/cents_kperm_stack.mat' ] 
save(filename,'data','-v7.3')
ls(filename)
   
%% segs_vox

data = zeros(NNODES,1,n_animal,NRESAMP) ;
for idx = 1:n_animal

    disp(idx)

    animal_str = all_data{idx} ;

    for jdx = 1:NRESAMP

        data(:,:,idx,jdx) = animal_str(jdx).dat.segs_vox ;

    end
end
filename = [ DD_INTERM '/segs_vox_kperm_stack.mat' ] 
save(filename,'data','-v7.3')
ls(filename)
   
%% sanity check

% sqr_entries = { 'dist_mat' 'leng_mat' 'node_adj' 'con_mat' 'con_mat_gn' } ;
%     
% % save data
% filename = [ DD_INTERM '/' sqr_entries{5} '_kperm_stack.mat' ] 
% ll = load(filename)
% ls(filename)
% 
% %% viz for sanity
% 
% % across animals
% imsc_overdim(ll.data(:,:,1:226),3,[],1)
% 
% % across reps
% imsc_overdim(squeeze(ll.data(:,:,25,:)),3,[],1)
% imsc_overdim(squeeze(ll.data(:,:,101,:)),3,[],1)


