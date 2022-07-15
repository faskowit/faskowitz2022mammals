% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

VIZIT = 0 ;

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

% filename = [ DD_INTERM '/cents_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
% ll = load(filename) ; 
% c_aCent = ll.data ;
% clear ll ;

% filename = [ DD_INTERM '/node_adj_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
% ll = load(filename) ; 
% c_aAdj = ll.data ;
% clear ll ;

%% setup

triunroll = @(x_) x_(logical(triu(ones(size(x_)),1))) ; 
randlevels = [ 0 ((5:5:100) ./ 100) ] ;
n_animal = size(ssheet,1) ;

%% init outputs

hub_resil = cell(n_animal,1) ;

hub_resilarea = nan(n_animal,1) ;
hub_resiljac = nan(n_animal,1) ;
hub_resilvar = nan(n_animal,1) ;

%% big loop

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_repanimal_hubresil_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
if ~isfile(filename)

%%
    
parfor idx = 1:n_animal

    disp(idx) 
    
    % data 
    c = c_aCon(:,:,idx) ;
%     cents = c_aCent(:,:,idx) ;
%     d = squareform(pdist(cents)) ;
    d = c_aLeng(:,:,idx) ;

    try

        % function [ hubresil, hubresilarea, resiljac, resilvar, orighub ] = ...
        % hub_resilience(mat,distmat,hubmethod,Nreps)
        [hr,hra,rj,rv] = hub_resilience(c,d,'mydef',250) ;
        
        hub_resil{idx} = hr(:) ;
        
        hub_resilarea(idx) = hra ;
        hub_resiljac(idx) = rj ;
        hub_resilvar(idx) = rv ; 

    catch
        warning('catch error') 

    end

end

else % if data is already generated, load it
   load(filename) 
end

%% saveit

filename = [ DD_PROC '/' OUTSTR '_repanimal_hubresil_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'hub_*','randlevels','-v7.3')


%% end the big threshold loop
end

%% some viz stuff

if VIZIT

thr=0.1 ;
filename = [ DD_PROC '/' OUTSTR '_repanimal_hubresil_thr' num2str(thr) '_.mat' ] ;
load(filename)
filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
ssheet = ll.newsheet ;

hubstack = zeros(200,length(hub_resil)) ;
for idx = 1:length(hub_resil)
   hubstack(:,idx) = hub_resil{idx}(~isnan(hub_resil{idx})) ;
end

nice_scatter(ssheet.log10_BrV_,mean(hubstack),200,grp2idx(ssheet.Order)) % all the same
nice_scatter(ssheet.log10_BrV_,median(hubstack),200,grp2idx(ssheet.Order)) 

nice_scatter(ssheet.log10_BrV_,std(hubstack),200,grp2idx(ssheet.Order)) % all the same

nice_scatter(ssheet.log10_BrV_,sum(hubstack==0),200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,sum(hubstack<0.1),200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,sum(hubstack>0.1),200,grp2idx(ssheet.Order)) 

for idx = 0:0.05:1
nice_scatter(ssheet.log10_BrV_,sum(hubstack>idx),200,grp2idx(ssheet.Order)) 
title(idx)
disp(idx)
waitforbuttonpress
end

end

%%

if VIZIT

    % mean across reps
    nice_scatter(ssheet.log10_BrV_,hub_resilarea,200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,hub_resiljac,200,grp2idx(ssheet.Order)) 
    nice_scatter(ssheet.log10_BrV_,hub_resilvar,200,grp2idx(ssheet.Order)) 
  
end


