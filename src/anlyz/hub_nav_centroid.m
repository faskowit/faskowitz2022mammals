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
excludiag = @(x_) x_(~logical(eye(size(x_,1)))) ;

% randlevels = [ 0 ((5:5:100) ./ 100) ] ;
n_animal = size(ssheet,1) ;

m = logical(triu(ones(NNODES),1)) ;

%% init outputs

nReps = 1000 ;

nav_eff = nan(n_animal,1) ;
nav_eff_hh = nan(n_animal,1) ;
nav_eff_nn = nan(n_animal,1) ;
dens = nan(n_animal,1) ;

hubscore = nan(n_animal,NNODES) ;

nav_cent_hh = nan(n_animal,1) ;
nav_cent_nn = nan(n_animal,1) ;

rep_naveff = nan(n_animal,nReps) ;
rep_naveff_hh = nan(n_animal,nReps) ;
rep_naveff_nn = nan(n_animal,nReps) ;

rep_navcent_hh = nan(n_animal,nReps) ;
rep_navcent_nn = nan(n_animal,nReps) ; 

%% big loop

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_hubnav_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
if ~isfile(filename)

%%
parfor idx = 1:n_animal

    disp(idx) 
    
    % data 
    c = c_aCon(:,:,idx) ;
    cents = c_aCent(:,:,idx) ;
    eud = squareform(pdist(cents)) ;
    d = c_aLeng(:,:,idx) ;
    mask = ~eye(NNODES) ;

%     conmask = c > 0 ;

    gc = get_components(c) ;
    if any(gc~=1)
        mm = mode(gc) ; 
        vv = mm(1) == gc ; 
        c = c(vv,vv) ;
        cents = cents(vv,:) ;
        eud = eud(vv,vv) ;
        d = d(vv,vv) ;
        mask = mask(vv,vv) ;
    end
   
    try

        dens(idx) = density_und(c) ;
        
        % Weighted 
        inv_wei = 1./c ;
        inv_wei(isinf(inv_wei)) = 0 ;
        
        % function [sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(L, D, max_hops)
        [~,plbin,~,~,pp] = navigation_wu(inv_wei,eud,100) ;
        spdist = distance_bin(c>0) ; 
        
        [hs,hh] = get_hub_score_wei_und(c,'mydef') ;
       
        hubscore(idx,:) = hs ;
        
        nav_ratio = spdist ./ plbin ; 
        nav_cent = nav_centrality(plbin,pp) ;

        nav_eff(idx) = mean(excludiag(nav_ratio)) ; 
        
        nav_cent_hh(idx) = mean(nav_cent(hh)) ;
        nav_cent_nn(idx) = mean(nav_cent(~hh)) ;

        nav_eff_hh(idx) = mean(excludiag(nav_ratio(hh,hh))) ;
        nav_eff_nn(idx) = mean(excludiag(nav_ratio(~hh,~hh))) ;
        
        %% rewire ish
        
        tic
        for ndx = 1:nReps
            
            % disp(ndx)
            
            [~,~,rewr_eud] = geombinsurr_partial(eud,eud,1,3,'quantiles') ;

            [~,rw_plbin,~,~,rwpp] = navigation_wu(inv_wei,rewr_eud,100) ;

            rw_nv = spdist ./ rw_plbin ; 
            rw_nc = nav_centrality(rw_plbin,rwpp) ;
            
            rep_naveff(idx,ndx) = mean(excludiag(rw_nv)) ; 
            rep_naveff_hh(idx,ndx) = mean(excludiag(rw_nv(hh,hh))) ;
            rep_naveff_nn(idx,ndx) = mean(excludiag(rw_nv(~hh,~hh))) ;
            
            rep_navcent_hh(idx,ndx) = mean(rw_nc(hh)) ;
            rep_navcent_nn(idx,ndx) = mean(rw_nc(~hh)) ;
            
        end
        disp(['nulls took: ' num2str(toc) ])
        
    catch
        warning('caught error') 
    end

%     try
% 
%         % unwighted 
%         un_wei = c>0 ;
% 
% 
%     catch
%         warning('caught error') 
%     end
%     
end

else % if data is already generated, load it
   load(filename) 
end

%% saveit

filename = [ DD_PROC '/' OUTSTR '_hubnav_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'*nav*','hubscore','dens','-v7.3')


%% end the big threshold loop
end

%% vizit

if VIZIT

thr=0;
filename = [ DD_PROC '/' OUTSTR '_hubnav_thr' num2str(thr) '_.mat' ] ;
load(filename)
filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
ssheet = ll.newsheet ;

nice_scatter(ssheet.log10_BrV_,nav_eff,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_eff_hh,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_eff_nn,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_eff_hh./(nav_eff_hh+nav_eff_nn),200,grp2idx(ssheet.Order)) 

nice_scatter(ssheet.log10_BrV_,nav_cent_hh+nav_cent_nn,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_cent_hh,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_cent_nn,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,nav_cent_hh./(nav_cent_hh+nav_cent_nn),200,grp2idx(ssheet.Order)) 


[~,~,resid] = regress(nav_eff,dens) ;
nice_scatter(ssheet.log10_BrV_,resid,200,grp2idx(ssheet.Order)) 


% nice_scatter(ssheet.log10_BrV_,,200,grp2idx(ssheet.Order)) 

end

%% permutations

nullcorr = nan(nReps,1) ;

for idx = 1:nReps
    
    nullcorr(idx) = corr(ssheet.log10_BrV_,rep_naveff(:,idx),'rows','complete','type','s') ;
    
end






