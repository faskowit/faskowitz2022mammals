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

% filename = [ DD_INTERM '/leng_mat_repani_stack_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
% ll = load(filename) ; 
% c_aLeng = ll.data ;
% clear ll ;

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

nparams = 100;

estK = nan(n_animal,1) ;
estPar = nan(n_animal,1) ;

%% big loop

% make sure the data isn't already generated
filename = [ DD_PROC '/' OUTSTR '_generative_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
if ~isfile(filename)

%%
parfor idx = 1:n_animal

    %%
    
    disp(idx) 
    
    % data 
    c = c_aCon(:,:,idx) ;
    cents = c_aCent(:,:,idx) ;
    eud = squareform(pdist(cents)) ;

    amst = graphminspantree(sparse(1./c),'method','kruskal');
    Cseed = full(amst + amst')>0;

%     % select 90% of MST
%     seedinit = find(amst) ;
%     seedind = randperm(length(seedinit),ceil(NNODES*.9)) ;
%     seeds = seedinit(seedind) ;
%     Cseed = zeros(size(c)) ;
%     Cseed(seeds) = 1 ;
%     Cseed = Cseed + Cseed' ; 
    
    % set model type
    modeltype = 'sptl';

    % set whether the model is based on powerlaw or exponentials
    modelvar = [{'powerlaw'},{'powerlaw'}];

    % choose some model parameters
    lowparam = -15 ;
    highparam = 0 ; 
    
    for ndx = 1:5
    
        tic
        disp(ndx) 
        
        srchparam = unifrnd(lowparam,highparam,nparams,1);

        % generate synthetic networks and energy for the neighbors model;
        [~,K] = evaluate_distonly_generative_model(Cseed,c,eud,modeltype,modelvar,srchparam);

        % sort K
        [sortedK,sortedKind] = sort(K,'ascend') ;

        % get the lower inds
        bottomPrct = .10 ;
        
        selectNum = floor(nparams*bottomPrct) ;
        expandNum = floor(nparams*bottomPrct*2) ;

        bottomInds = sortedKind(1:selectNum) ;
                
        % now sample with bias towards lowest energies
        bti2 = randsample(expandNum,selectNum,true,(1-sortedK(1:expandNum))) ;
        bottomInds2 = sortedKind(bti2) ;
        
        comboBottom = unique([ bottomInds ; bottomInds2 ])
        
        lowparam = min(srchparam(comboBottom)) 
        highparam = max(srchparam(comboBottom)) 
        toc
        
    end
    
    
    % estK(idx) = quantile(K,0.05) ;
    
    finalinds = sortedKind(1:bottomInds) ;
       
    % take the median of the bottom inds
    med = median(K(finalinds),1) ;
    [~,medInds] = min(abs(med-K(finalinds))) ;
    
    estK(idx) = K(finalinds(medInds(1))) ;
    estPar(idx) = srchparam(finalinds(medInds(1))) ; 
    
%     scatter(srchparam,K,100,K,'filled') ; 
%     set(gca,...
%     'ylim',[0,1],...
%     'clim',[0,1]);
%     colormap(turbo);
%     xlabel('geometric parameter, \eta');
%     ylabel('edge length');
      
end

else % if data is already generated, load it
   load(filename) 
end

%% saveit

filename = [ DD_PROC '/' OUTSTR '_generative_thr' num2str(THRDENS(tDx)) '_.mat' ] ;
save(filename,'estK','estPar','-v7.3')


%% end the big threshold loop
end

%%

if VIZIT

thr=0.1;
filename = [ DD_PROC '/' OUTSTR '_generative_thr' num2str(thr) '_.mat' ] ;
load(filename)
filename = [ DD_INTERM '/con_mat_gn_repani_stack_thr' num2str(thr) '_.mat' ] ;
ll = load(filename,'newsheet') ; 
ssheet = ll.newsheet ;

nice_scatter(ssheet.log10_BrV_,estPar,200,grp2idx(ssheet.Order)) 
nice_scatter(ssheet.log10_BrV_,estK,200,grp2idx(ssheet.Order)) 

nice_scatter(ssheet.log10_BrV_,ssheet.SNR,200,grp2idx(ssheet.Order)) 

[~,~,resid] = regress(estPar,ssheet.SNR) ;
nice_scatter(ssheet.log10_BrV_,resid,200,grp2idx(ssheet.Order)) 

end
