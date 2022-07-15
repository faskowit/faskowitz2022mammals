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
% n_animal = length(animal_names) ;

VIZIT = 0 

%% big loop

for tDx = 1:length(THRDENS)
%%

filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
ll = load(filename) ; 
c_aCon = ll.data ;
ssheet = ll.newsheet ;
n_animal = size(ssheet,1) ;
clear ll ;

filename = [ DD_INTERM '/leng_mat_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
ll = load(filename) ; 
c_aLeng = ll.data ;
clear ll ;

%% init

randlevels = [ 0 ((5:5:100) ./ 100) ] ;

geo_dist = cell(n_animal,1) ;
rand_dist = cell(n_animal,1) ;
randmio_dist = cell(n_animal,1) ;

%% distances yo

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_repanimal_rewrdist_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
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

        nreps = 100 ;
        % nreps = 10 ;
        
%         function [geo_netdist,rand_netdist] = ...
%             net_rewire_dist(mat,distmat,Nreps,prop_rand,nbins)
        gd = netpd_rewire_dist(c,d,'gs',nreps,randlevels) ;
        rd = netpd_rewire_dist(c,d,'sw',nreps,randlevels) ;
        ra = netpd_rewire_dist(c,d,'rm',nreps,randlevels) ;
        
        geo_dist{idx} = gd ;
        rand_dist{idx} = rd ;
        randmio_dist{idx} = ra ;

    catch
        warning('catch error') 
        
        geo_dist{idx} = nan ;
        rand_dist{idx} = nan ;
        randmio_dist{idx} = nan ;
    end

end

else % if data is already generated, load it
   load(filename) 
end

% save it

filename = [ DD_PROC '/' OUTSTR '_repanimal_rewrdist_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'*_dist','randlevels','-v7.3')

end

%%

if VIZIT

%     % 100% rewire
%     nice_scatter(ssheet.log10_BrV_,...
%         cellfun(@(x_)x_(end,1),geo_dist ),200,grp2idx(ssheet.Order)) 
%     ylim([0 0.5]) ; xlim([1 6])
%     nice_scatter(ssheet.log10_BrV_,...
%         cellfun(@(x_)x_(end,1),rand_dist ),200,grp2idx(ssheet.Order)) 
%     ylim([0 0.5]) ; xlim([1 6])
%     nice_scatter(ssheet.log10_BrV_,...
%         cellfun(@(x_)x_(end,1),randmio_dist ),200,grp2idx(ssheet.Order)) 
%     ylim([0 0.5]) ; xlim([1 6])
% 
%     nice_scatter(...
%         cellfun(@(x_)x_(11,1),rand_dist ),...
%         cellfun(@(x_)x_(11,1),geo_dist ),200,grp2idx(ssheet.Order)) 
%     hold on 
%     refline(1,0)
%     hold off
%     
%     % 5% rewire
%     nice_scatter(ssheet.log10_GM_,...
%         cellfun(@(x_)x_(2,1),geo_dist ),200,grp2idx(ssheet.Order))   
%     nice_scatter(ssheet.log10_GM_,...
%         cellfun(@(x_)x_(2,1),rand_dist ),200,grp2idx(ssheet.Order)) 
%     
%     % 50% rewire
%     nice_scatter(ssheet.log10_GM_,...
%         cellfun(@(x_)x_(11,1),geo_dist ),200,grp2idx(ssheet.Order))   
%     nice_scatter(ssheet.log10_GM_,...
%         cellfun(@(x_)x_(11,1),rand_dist ),200,grp2idx(ssheet.Order)) 
    
    
    subplot(1,3,1) 
    nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),geo_dist ),...
        200,grp2idx(ssheet.Order)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 
   
    subplot(1,3,2) 
    nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),rand_dist ),...
        200,grp2idx(ssheet.Order)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 

    subplot(1,3,3) 
    nice_scatter(ssheet.log10_BrV_,...
        cellfun(@(x_)trapz(randlevels,x_(:,1)./max(x_(:,1))),randmio_dist ),...
        200,grp2idx(ssheet.Order)) 
    xlim([ 1.5 6.5 ]) ; ylim( [ 0.5 1 ] ) 


    for idx = 2:21

        nice_scatter(...
            cellfun(@(x_)x_(idx,1),rand_dist ),...
            cellfun(@(x_)x_(idx,1),geo_dist ),200,grp2idx(ssheet.Order)) 
%         nice_scatter(...
%             cellfun(@(x_)x_(idx,1),rand_dist ),...
%             cellfun(@(x_)x_(idx,1),geo_dist ),200,ones(n_animal,1).*idx) 
        title([ 'randomize: ' num2str(randlevels(idx)) ])
        xlim([0 0.5])
        ylim([0 0.5])
        hold on 
        refline(1,0)
%         hold off 

        xlabel('random rewire')
        ylabel('geo rewire')
        
%         waitforbuttonpress

    end

    
end

%%

aa = grp2idx(ssheet.Order) ;
n_ani_classes = length(unique(grp2idx(ssheet.Order))) ;
ani_colors = parula(n_ani_classes) ;
ani_names = unique(ssheet.Order) ;

for idx = 1:n_ani_classes

    disp(idx)
    
    tmpcell = geo_dist(aa==idx) ;
    
    subplot(3,4,idx) 
    for jdx = 1:size(tmpcell,1) 

        col = ani_colors(idx,:) ;
        
        hold on
        % function [] = plot_smokey(xvals,yvals,ystd,color1,color2,plotpattern)
        plot_smokey(randlevels,tmpcell{jdx}(:,1),tmpcell{jdx}(:,2),...
            col,col,[1 1])
        ylim([ 0 0.4])        
    end
            title(ani_names{idx})
    
end

hold off

