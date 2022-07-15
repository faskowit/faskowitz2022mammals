% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

filename = [ DD_INTERM '/con_mat_gn_kperm_stack.mat' ] ;
ll = load(filename) ;
allCon = ll.data ; clear ll ;

%%

n_animals = size(allCon,3) ;
allComms = cell(n_animals,1) ;

for idx = 1:n_animals
    
    animalCon = squeeze(allCon(:,:,idx,:)) ;
        
    tempRes = struct() ;

    parfor jdx = 1:NRESAMP
        
        tmpCon = animalCon(:,:,jdx) ;
                
        %       K:      network data
        %       sam1:   samples communities first loop
        %       sam2:   samples communities second loop
        %       maxC:   max num communities
        %       r1:     initial gamma low
        %       r2:     initial gamma high
        %       tau:    tau parameter for consensus modules at the end 
        
        %       ciu:    consensus community partition
        %       Aall:   agreement matrix
        %       anull:  analytic null
        %       A:      agree - null
        %       ciall:  community assignments
        
        samp1 = 500;
        samp2 = 2500;
        N = 100;
        r1 = -1;
        r2 = 3;

        [ciu,~,~,A] = get_SCmodules_MRCClite(tmpCon,samp1,samp2,N,r1,r2) ;
        
        tempRes(jdx).ciu = ciu ;
        tempRes(jdx).A = A ;
        
    end
    
    allComms{idx} = tempRes ; 
    
    disp([ 'FINISHED ANIMAL: ' num2str(idx) ])
    
end

%% reformat struct

data = struct() ;

for idx = 1:n_animals
   
    disp(idx) 
    
    data(idx).ciu = zeros(NNODES,NRESAMP) ;
    data(idx).Amat = zeros(NNODES,NNODES,NRESAMP) ;
    
    for jdx = 1:NRESAMP
        
        data(idx).ciu(:,jdx) = allComms{idx}(jdx).ciu ;
        data(idx).Amat(:,:,jdx) = allComms{idx}(jdx).A ;
        
    end
end

%% save eet

filename = [ DD_INTERM '/' OUTSTR '_cons_comms.mat' ] ;
save(filename,'data','-v7.3')
ls(filename)




