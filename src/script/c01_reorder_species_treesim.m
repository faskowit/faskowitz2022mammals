% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

filename = [ DD_PROC '/MamiNum_corrected_proc_2.csv' ] ;
ssheet = readtable(filename,'Delimiter','\t') ;

%%

rt = readtable([ DD_RAW '/vertLife/proc/vertLife_distmat.txt'],...
                 'ReadRowNames',1,'ReadVariableNames',1) ;
             
spMat = table2array(rt) ; 
             
rt_rownames = rt.Properties.RowNames ;  
rt_rownames = lower(strrep(rt_rownames,' ','_')) ;

sh_rownames = ssheet.sci_name ;
sh_rownames = lower(strrep(sh_rownames,' ','_')) ;

%% now match it up
           
[~,an2spmatInd]= ismember(sh_rownames,rt_rownames) ;
an2spmatIndNoZ = an2spmatInd(an2spmatInd>0) ;

%% make new spMat 

speciesSim = nan(length(sh_rownames)) ;
speciesSim(an2spmatInd>0,an2spmatInd>0) = spMat(an2spmatIndNoZ,an2spmatIndNoZ) ;

%%

filename = [ DD_PROC '/species_treesim.mat'];
save(filename,'speciesSim','an2spmatInd','-v7.3')


