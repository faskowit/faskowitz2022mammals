% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% read in proc spreadsheet

% filename = [ DD_PROC '/MamiNum_corrected_proc.txt' ] ;
% ssheet = readtable(filename,'Delimiter','\t') ;
% 
% animal_names = regexprep(ssheet.file_prefix,'_FF3_allfibers_kperm','') ;
% % n_animal = length(animal_names) ;

VIZIT = 0 

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

% filename = [ DD_INTERM '/node_adj_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
% ll = load(filename) ; 
% c_aAdj = ll.data ;
% clear ll ;

%% setup

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 
randlevels = [ 0 ((5:5:100) ./ 100) ] ;
n_animal = size(ssheet,1) ;

%% init outputs

hub_prob = cell(n_animal,1) ;
hub_frag = cell(n_animal,1) ;
hub_dist = nan(n_animal,1) ;
hub_sus = nan(n_animal,length(randlevels)) ;
hub_entropy = nan(n_animal,length(randlevels)) ;
hub_ent_trap = nan(n_animal,1) ;

%% big loop

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_fragilityMyDef2_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
if ~isfile(filename)

%%
parfor idx = 1:n_animal

    disp(idx) 
    
    % data 
    c = c_aCon(:,:,idx) ;
    cents = c_aCent(:,:,idx) ;
    eud = squareform(pdist(cents)) ;
    d = c_aLeng(:,:,idx) ;
    
    % clean up mat
    gc = get_components(c) ;
    gcmode = mode(gc) ;
    if any(gc~=gcmode(1))
        % make a lil smaller 
        eud = eud(gc==gcmode(1),gc==gcmode(1)) ;
        c = c(gc==gcmode(1),gc==gcmode(1)) ;
        d = d(gc==gcmode(1),gc==gcmode(1)) ; 
    end    

    try

        [orighub,hubprob,hubarea,hubsus] = ...
            hub_sweep(c,d,'mydef2',250,randlevels,eud) ;
         
        % function [hubfrag] = get_hub_fragility(hubs,hubprob,proprand) 
        hubfrag = get_hub_fragility(orighub,hubprob,randlevels) ;
        
        % record fragility
        hub_frag{idx} = hubfrag(:) ;
        % susceptibility
        hub_sus(idx,:) = hubsus ;
        % entropy on hub probabilities
        hubp_norm = hubprob ./ sum(hubprob) ; % normalize
        ee = -sum(hubp_norm .* log2(hubp_norm),1,'omitnan') ./ log2(size(hubp_norm,1)) ;
        hub_entropy(idx,:) = ee ;
        hub_ent_trap(idx) = trapz(randlevels,ee) ; % entropy under curve
        
        hub_dist(idx) = trapz(randlevels,hubarea) ; % distance under curve

        hub_prob{idx} = hubprob ;
        
    catch
        warning('caught error') 

    end

end

else % if data is already generated, load it
   load(filename) 
end

%% saveit

filename = [ DD_PROC '/' OUTSTR '_fragilityMyDef2_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'hub_*','randlevels','-v7.3')


%% end the big threshold loop
end

%% some viz stuff

thr=0.1
filename = [ DD_PROC '/' OUTSTR '_fragilityMyDef2_thr' num2str(thr) '_.mat' ] ;
load(filename)
filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
ssheet = ll.newsheet ;
n_animal = size(ssheet,1) ;

% make hubfrag2
hubfragAUC = cell(n_animal,1) ;
for idx = 1:n_animal
    [~,hubfragAUC{idx}] = get_hub_fragility(hub_prob{idx}(:,1),hub_prob{idx},randlevels) ;
end


hfstack = zeros(30,length(hub_frag)) ;
for idx = 1:length(hub_frag)
   hfstack(:,idx) = hub_frag{idx}(~isnan(hub_frag{idx})) ;
end

hfstack2 = zeros(30,length(hubfragAUC)) ;
for idx = 1:length(hubfragAUC)
   hfstack2(:,idx) = hubfragAUC{idx}(~isnan(hubfragAUC{idx})) ;
end

if VIZIT

    
    nice_scatter(ssheet.log10_BrV_,mean(hfstack),200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,mean(hfstack2),200,grp2idx(ssheet.Order)) 

    
    nice_scatter(ssheet.log10_BrV_,median(hfstack),200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,std(hfstack),200,grp2idx(ssheet.Order)) 

    nice_scatter(ssheet.log10_BrV_,sum(hfstack<0.1),200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,sum(hfstack>0.1),200,grp2idx(ssheet.Order)) 

end

%%

if VIZIT

    % mean across reps
    nice_scatter(ssheet.log10_BrV_,hub_dist,200,grp2idx(ssheet.Order))  %#ok<UNRCH>

    nice_scatter(ssheet.log10_BrV_,mean(hub_sus,2),200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,median(hub_sus,2),200,grp2idx(ssheet.Order)) 
    
    nice_scatter(ssheet.log10_BrV_,median(hub_entropy,2),200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,hub_ent_trap,200,grp2idx(ssheet.Order)) 
end






