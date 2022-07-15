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
n_animal = size(ssheet,1) ;
thr_vals = [ 0 0.05 0.1 0.15 ] ; 

%% load data

filename = [ DD_INTERM '/con_mat_gn_kperm_stack.mat' ] ;
ll = load(filename) ;
ls(filename)

%% add path to netpd

addpath(genpath('/home/spornslab/joshstuff/git_pull/netpd_matlab'))

%% first check if mat is fully connected
% information not used

% notconnected = zeros(n_animal,NRESAMP) ; 
% 
% for idx = 1:n_animal
%     
%     disp(idx)
%     
%     aCon = squeeze(ll.data(:,:,idx,:)) ; 
% 
%     for jdx = 1:NRESAMP
% 
%         c = aCon(:,:,jdx) ;
%         
%         gc = get_components(c) ;
%         gcmode = mode(gc) ;
%         if any(gc~=gcmode(1))
%             disp([ num2str(jdx) ' FEFFEFEFFEF' ])
%             notconnected(idx,jdx) = sum(gc~=gcmode(1)) ;
%         end
% 
%     end 
% end

%% compare networks across thresholds (LONG COMPUTE TIME)


netdivergences = cell(n_animal,length(thr_vals)) ;

for tdx = 1:length(thr_vals)

    currthr = thr_vals(tdx) ;
    
for idx = 1:n_animal
    
    divs = nan(NRESAMP) ;
   
    aCon = squeeze(ll.data(:,:,idx,:)) ; 
    
    parfor jdx = 1:NRESAMP
        
        disp([ 'thr: ' num2str(currthr) ' : ' num2str(idx) ' - ' num2str(jdx) ])
        
        c_j = aCon(:,:,jdx) ;
        if currthr > 0
            try
            	c_j = c_j .* full(fcn_mst_plus(c_j,currthr,'inv')) ;
            catch
                continue
            end
        end
        
        D1 = distance_wei_floyd(c_j,'inv') ;
        
        for kdx = 1:NRESAMP
            if jdx>=kdx ; continue ; end
                        
            c_k = aCon(:,:,kdx) ;
            if currthr > 0
                try
                    c_k = c_k .* full(fcn_mst_plus(c_k,currthr,'inv')) ;
                catch
                    continue
                end
            end
            
            D2 = distance_wei_floyd(c_k,'inv') ;
            
            % get commoin edge bins for distances
            ebins = netpd_edgevalbins(D1,D2,'quantiles',25) ;

            % two portraits, with common bins
            B1 = netpd_portrait_wei(D1,ebins,'alreadydistance') ;
            B2 = netpd_portrait_wei(D2,ebins,'alreadydistance') ;

            divs(jdx,kdx) = netpd_divergence(B1,B2) ;

        end
    end

    netdivergences{idx,tdx} = divs ;

end

end

%% throw away any major outliers within each mammal (i.e. bad random parcs)
% we do this, so that the bad network distances doesn't influence the
% selection of the centroid parition for each animal

% % viz it
% for idx = 1:n_animal
%    
%     for ndx = 1:4
%         subplot(1,4,ndx)
%         imagesc(netdivergences{idx,ndx},[0 0.15])
%     end
%     title([ num2str(idx) ' ' ssheet.AnimalName(idx)])
%     waitforbuttonpress
% end

netdivergences_cln = cell(n_animal,length(thr_vals)) ;

for tdx = 1:length(thr_vals)
    
for idx = 1:n_animal
   
    disp([ num2str(tdx) ' : ' num2str(idx) ])
        
    mat = netdivergences{idx,tdx} ;

    % make lower triangle 0
    mat(~~tril(ones(size(mat,1)),-1)) = 0 ;
    mat = mat + mat' ;
    
    dd = mean(mat,2,'omitnan') ;
    [~,~,oU] = isoutlier(dd) ; % only need upper outlier for now..
    new_mat = mat(dd<oU,dd<oU) ;
    
    netdivergences_cln{idx,tdx} = new_mat ;
end

end

%% find the centroid network, for each animals, at each threshold

leastdist = zeros(n_animal,length(thr_vals)) ;
leastdist2nd = zeros(n_animal,length(thr_vals)) ;
leastdist3rd = zeros(n_animal,length(thr_vals)) ;
verydist = zeros(n_animal,length(thr_vals)) ;

netdiv_distances = nan(NRESAMP,n_animal,length(thr_vals)) ;

triunroll = @(x_) x_(logical(triu(ones(size(x_,1)),1))) ;

for tdx = 1:length(thr_vals)

for idx = 1:n_animal
   
    disp([ num2str(tdx) ' : ' num2str(idx) ])
    
    currmat = netdivergences_cln{idx,tdx} ;
    
    if isempty(currmat)
       
        netdiv_distances(:,idx,tdx) = nan(NRESAMP,1) ;
        leastdist(idx,tdx) = nan ;
        leastdist2nd(idx,tdx) = nan ;
        leastdist3rd(idx,tdx) = nan ;
        verydist(idx,tdx) = nan ;
        continue
    end
    
    mat = currmat .* 1 ;
    % make lower triangle 0
    mat(~~tril(ones(size(mat,1)),-1)) = 0 ;
    mat = mat + mat' ;
    
    % distances
    dd = mean(mat,2,'omitnan') ;
    [dds,ss] = sort(dd,'ascend') ;
    
    netdiv_distances(1:length(dds),idx,tdx) = dds ;
    leastdist(idx,tdx) = ss(1) ;
    leastdist2nd(idx,tdx) = ss(2) ;
    leastdist3rd(idx,tdx) = ss(3) ;
    try
        verydist(idx,tdx) = ss(50) ;
    catch
        verydist(idx,tdx) = ss(end) ;
    end
end

end

%% viz animals distances all together

% imagesc(netdiv_distances(:,:,1),[0 .25])
% imagesc(netdiv_distances(:,:,2),[0 .25])
% imagesc(netdiv_distances(:,:,3),[0 .25])
% imagesc(netdiv_distances(:,:,4),[0 .25])
% can see that some have way higher distances genereally (columns) than
% compared to other animals... probably should mark these animals for
% exclusion because something likely wrong with the scan

o_vec_1 = zeros(n_animal,1) ;

for idx = 1:4

    m = mean(netdiv_distances(:,:,1),'omitnan') ;
    [~,~,oU] = isoutlier(m(:)) ;
    o_vec_1 = o_vec_1 | m(:) > oU ;
end

% viz it now
% imagesc(netdiv_distances(:,~o_vec_1,1),[0 .25])

%% least dist mats

leastdist_con_mat_gn = struct() ;
leastdist_con_mat_gn_2nd = struct() ;
leastdist_con_mat_gn_3rd = struct() ;
verydist_con_mat_gn  = struct() ;

for tdx = 1:length(thr_vals)

leastdist_con_mat_gn(tdx).mat = zeros(NNODES,NNODES,n_animal) ;
leastdist_con_mat_gn_2nd(tdx).mat = zeros(NNODES,NNODES,n_animal) ;
leastdist_con_mat_gn_3rd(tdx).mat = zeros(NNODES,NNODES,n_animal) ;
verydist_con_mat_gn(tdx).mat = zeros(NNODES,NNODES,n_animal) ;

currthr = thr_vals(tdx) ;    

for idx = 1:n_animal
    
    disp(idx)
    
    %% 1st
    
    ldind = leastdist(idx,tdx) ;
    
    if isnan(ldind)
        dat = nan(NNODES) ;
    else
        dat = ll.data(:,:,idx,ldind) ;
    end
    
    try
        dat = dat .* full(fcn_mst_plus(dat,currthr,'inv')) ;
    catch
        warning('NOT ENOUGH EDGES')
        dat = nan(NNODES) ; 
    end
        
    leastdist_con_mat_gn(tdx).mat(:,:,idx) = dat ;
    
    %% 2nd
    
    ldind = leastdist2nd(idx,tdx) ;
    
    if isnan(ldind)
        dat = nan(NNODES) ;
    else
        dat = ll.data(:,:,idx,ldind) ;
    end
    
    try
        dat = dat .* full(fcn_mst_plus(dat,currthr,'inv')) ;
    catch
        warning('NOT ENOUGH EDGES')
        dat = nan(NNODES) ; 
    end
        
    leastdist_con_mat_gn_2nd(tdx).mat(:,:,idx) = dat ;
       
    %% 3rd
    
    ldind = leastdist3rd(idx,tdx) ;
    
    if isnan(ldind)
        dat = nan(NNODES) ;
    else
        dat = ll.data(:,:,idx,ldind) ;
    end
    
    try
        dat = dat .* full(fcn_mst_plus(dat,currthr,'inv')) ;
    catch
        warning('NOT ENOUGH EDGES')
        dat = nan(NNODES) ; 
    end
        
    leastdist_con_mat_gn_3rd(tdx).mat(:,:,idx) = dat ; 
    
    %% very
    
    ldind = verydist(idx,tdx) ;
    
    if isnan(ldind)
        dat = nan(NNODES) ;
    else
        dat = ll.data(:,:,idx,ldind) ;
    end
    
    try
        dat = dat .* full(fcn_mst_plus(dat,currthr,'inv')) ;
    catch
        warning('NOT ENOUGH EDGES')
        dat = nan(NNODES) ; 
    end
        
    verydist_con_mat_gn(tdx).mat(:,:,idx) = dat ; 
end

end

%% pick the most central
% 
% sort_inds = zeros(n_animal,NRESAMP) ;
% rep_mats_leastdist = zeros(NNODES,NNODES,n_animal) ;
% rep_mats_mostdist = zeros(NNODES,NNODES,n_animal) ;
% 
% for idx = 1:n_animal
%     
%     curr_mat = netdivergences{idx} ;
%    
%     curr_mat = curr_mat + curr_mat' ;
%     
%     [~,sIdx] = sort(sum(curr_mat),'ascend') ;
%     
%     sort_inds(idx,:) = sIdx ;
%     
%     rep_mats_leastdist(:,:,idx) = ll.data(:,:,idx,sort_inds(idx,1)) ;
%     rep_mats_mostdist(:,:,idx) = ll.data(:,:,idx,sort_inds(idx,100)) ;
% 
% end

%% make an animal divergence matrix (takes a minute or two)

animal_div_mat = nan(n_animal,n_animal,length(thr_vals)) ;

for tdx = 1:length(thr_vals)
    
    aCon = leastdist_con_mat_gn(tdx).mat ;
    
parfor idx = 1:n_animal
    
    i_dist = distance_wei_floyd(aCon(:,:,idx),'inv') ;
    
    for jdx = 1:n_animal
        if idx>=jdx ; continue ; end
        
        disp([ 'thr: ' num2str(tdx) ' : ' num2str(idx) ' - ' num2str(jdx) ])

        j_dist = distance_wei_floyd(aCon(:,:,jdx),'inv') ;
        
        if ( sum(i_dist,'all','omitnan') == 0 ) || ...
                (sum(j_dist,'all','omitnan') == 0)
            continue
        end
            
        % get commoin edge bins for distances
        ebins = netpd_edgevalbins(i_dist,j_dist,'quantiles',25) ;

        % two portraits, with common bins
        B1 = netpd_portrait_wei(i_dist,ebins,'alreadydistance') ;
        B2 = netpd_portrait_wei(j_dist,ebins,'alreadydistance') ;
        
        animal_div_mat(idx,jdx,tdx) = netpd_divergence(B1,B2) ;
        
    end
end

end

%% make the animal div mats for the 2nd and 3rd as well

animal_div_mat_2nd = nan(n_animal,n_animal,length(thr_vals)) ;
animal_div_mat_3rd = nan(n_animal,n_animal,length(thr_vals)) ;

for tdx = 1:length(thr_vals)
    
    aCon = leastdist_con_mat_gn_2nd(tdx).mat ;
    
parfor idx = 1:n_animal
    
    i_dist = distance_wei_floyd(aCon(:,:,idx),'inv') ;
    
    for jdx = 1:n_animal
        if idx>=jdx ; continue ; end
        
        disp([ 'thr: ' num2str(tdx) ' : ' num2str(idx) ' - ' num2str(jdx) ])

        j_dist = distance_wei_floyd(aCon(:,:,jdx),'inv') ;
        
        if ( sum(i_dist,'all','omitnan') == 0 ) || ...
                (sum(j_dist,'all','omitnan') == 0)
            continue
        end
            
        % get commoin edge bins for distances
        ebins = netpd_edgevalbins(i_dist,j_dist,'quantiles',25) ;

        % two portraits, with common bins
        B1 = netpd_portrait_wei(i_dist,ebins,'alreadydistance') ;
        B2 = netpd_portrait_wei(j_dist,ebins,'alreadydistance') ;
        
        animal_div_mat_2nd(idx,jdx,tdx) = netpd_divergence(B1,B2) ;
        
    end
end

end

for tdx = 1:length(thr_vals)
    
    aCon = leastdist_con_mat_gn_3rd(tdx).mat ;
    
parfor idx = 1:n_animal
    
    i_dist = distance_wei_floyd(aCon(:,:,idx),'inv') ;
    
    for jdx = 1:n_animal
        if idx>=jdx ; continue ; end
        
        disp([ 'thr: ' num2str(tdx) ' : ' num2str(idx) ' - ' num2str(jdx) ])

        j_dist = distance_wei_floyd(aCon(:,:,jdx),'inv') ;
        
        if ( sum(i_dist,'all','omitnan') == 0 ) || ...
                (sum(j_dist,'all','omitnan') == 0)
            continue
        end
            
        % get commoin edge bins for distances
        ebins = netpd_edgevalbins(i_dist,j_dist,'quantiles',25) ;

        % two portraits, with common bins
        B1 = netpd_portrait_wei(i_dist,ebins,'alreadydistance') ;
        B2 = netpd_portrait_wei(j_dist,ebins,'alreadydistance') ;
        
        animal_div_mat_3rd(idx,jdx,tdx) = netpd_divergence(B1,B2) ;
        
    end
end

end

%% very dist

animal_div_mat_very = nan(n_animal,n_animal,length(thr_vals)) ;

for tdx = 1:length(thr_vals)
    
    aCon = verydist_con_mat_gn(tdx).mat ;
    
parfor idx = 1:n_animal
    
    i_dist = distance_wei_floyd(aCon(:,:,idx),'inv') ;
    
    for jdx = 1:n_animal
        if idx>=jdx ; continue ; end
        
        disp([ 'thr: ' num2str(tdx) ' : ' num2str(idx) ' - ' num2str(jdx) ])

        j_dist = distance_wei_floyd(aCon(:,:,jdx),'inv') ;
        
        if ( sum(i_dist,'all','omitnan') == 0 ) || ...
                (sum(j_dist,'all','omitnan') == 0)
            continue
        end
            
        % get commoin edge bins for distances
        ebins = netpd_edgevalbins(i_dist,j_dist,'quantiles',25) ;

        % two portraits, with common bins
        B1 = netpd_portrait_wei(i_dist,ebins,'alreadydistance') ;
        B2 = netpd_portrait_wei(j_dist,ebins,'alreadydistance') ;
        
        animal_div_mat_very(idx,jdx,tdx) = netpd_divergence(B1,B2) ;
        
    end
end

end

%% now find outliers based on between-animal distances

% ca_ani = grp2idx(ssheet.Order) ;
% dat = animal_div_mat(:,:,1) ; 
% dat(~~tril(ones(size(dat,1)),-1)) = 0 ;
% dat = dat + dat' ; 
% dat = dat(~o_vec_1,~o_vec_1) ;
% ca_ani = ca_ani(~o_vec_1) ;
% 
% rr = 1:length(dat) % reorder_mod(dat,ca_ani) ;
% imsc_grid_comm(dat(rr,rr),ca_ani(rr),[],[],1)
% caxis([0 0.5])
% % clearly some weird animals here... like some animals that are distant
% % from all other animals, but especially within their own species family

outliers_thr = zeros(n_animal,length(thr_vals)) ;
outliers_thr_nonan = zeros(n_animal,length(thr_vals)) ;

for tdx = 1:4

    dat = animal_div_mat(:,:,tdx) ; 
    dat(~~tril(ones(size(dat,1)),-1)) = 0 ;
    dat = dat + dat' ;
    dat(1:size(dat,1)+1:end) = 0 ;

    % exclude the nanrows that pull the mean down
    nan_vec = sum(dat,2,'omitnan')==0 ;
    
    ca_ani = grp2idx(ssheet.Order) ;

    mm = mean(dat(~nan_vec,~nan_vec),2,'omitnan') ;
    mdz = module_degree_zscore(dat(~nan_vec,~nan_vec),ca_ani(~nan_vec)) ;

    [~,~,oU_1] = isoutlier(mm(:)) ;
%     [~,~,oU_2] = isoutlier(mdz(:)) ;
    oU_2 = 3 ; % positive standard deviations away from mean

    oo1 = zeros(n_animal,1) ;
    oo2 = zeros(n_animal,1) ;
    
    oo1(~nan_vec) = mm(:) > oU_1 ;
    oo2(~nan_vec) = mdz(:) > oU_2 ;
    
    outliers_thr(:,tdx) = o_vec_1 | oo1 | oo2 ;
    outliers_thr_nonan(:,tdx) =  outliers_thr(:,tdx) | nan_vec ;
end

%% visualize the divergence mats

% for tdx = 1:4
%     
%     outl = outliers_thr(:,tdx) ;
%     
%     dat = animal_div_mat(:,:,tdx) ; 
%     dat(~~tril(ones(size(dat,1)),-1)) = 0 ;
%     dat = dat + dat' ;
%     dat(1:size(dat,1)+1:end) = 0 ;
% 
%     % exclude the nanrows that pull the mean down
%     nan_vec = sum(dat,2,'omitnan')==0 ;
% 
%     keep = ~nan_vec & ~outl ; 
%     
%     ca_ani = grp2idx(ssheet.Order) ;
% 
%     imsc_grid_comm(dat(keep,keep),ca_ani(keep))
%     caxis([0 0.5])
%     
%     waitforbuttonpress
%     
% end

%% save it

filename = [ DD_INTERM '/' OUTSTR '_divergences_n_outlierdetect.mat' ] ;
save(filename,'netdivergences',...
              'netdivergences_cln',...
              'netdiv_distances',...
              'animal_div_mat*',...
              'leastdist*',...
              'verydist*',...
              'outliers_thr','outliers_thr_nonan','-v7.3') ;

filename = [ DD_INTERM '/' OUTSTR '_just_cent_ani_dat.mat' ] ;
save(filename,'leastdist*',...
              'verydist*',...
              'outliers_thr','outliers_thr_nonan','-v7.3') ;
