% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% inital read

% opts = detectImportOptions('./data/raw/MamiNum.xlsx');
orig_xls = readtable('./data/raw/MamiNum.xlsx') ; 

animal_names = orig_xls.AnimalName ; 
% cleanup
animal_names = regexprep(animal_names,'\(.*\)','') ;
animal_names = regexprep(animal_names,'_','') ;
animal_names = regexprep(animal_names,'\s','') ;
animal_names = lower(animal_names) ;

%% get list of files

all_mats = dir([ DD_RAW  '/xkperm/*kperm1.mat'  ]) ;

%% check 

warn = 0 ;
for idx = 1:length(all_mats)
    
    curr = all_mats(idx) ;
    
    a_name = regexprep(curr.name,'_FF3.*','') ;
    a_name = regexprep(a_name,'\(.*\)','') ;
    a_name = regexprep(a_name,'_','') ;
    a_name = regexprep(a_name,'\s','') ;
    a_name = lower(a_name) ;

    ff = strcmp(a_name,animal_names) ;
    
    if sum(ff) ~= 1
       warning(['fwffwf: ' a_name ]) 
       warn = warn +1 ;
       % waitforbuttonpress
    end
    
end

% at this point, add new column of corrected names and save a new csv

%% read revised

corr_xls = readtable('./data/raw/MamiNum_corrected.xlsx') ; 

animal_names_rev = corr_xls.CorrectedAnimalName ; 
% cleanup
animal_names_rev = regexprep(animal_names_rev,'\(.*\)','') ;
animal_names_rev = regexprep(animal_names_rev,'_','') ;
animal_names_rev = regexprep(animal_names_rev,'\s','') ;
animal_names_rev = lower(animal_names_rev) ;

%% figure out which of the 226 is not here ...

findInds = zeros(length(all_mats),1) ;
fn_pre = cell(length(all_mats),1) ;

for idx = 1:length(all_mats)
    
    curr = all_mats(idx) ;
    
    a_name = regexprep(curr.name,'_FF3.*','') ;
    a_name = regexprep(a_name,'\(.*\)','') ;
    a_name = regexprep(a_name,'_','') ;
    a_name = regexprep(a_name,'\s','') ;
    a_name = lower(a_name) ;

    findInds(idx) = find(strcmp(a_name,animal_names_rev)) ;
    
    % also get the prefixes of the filenames
    fn_pre{idx} = regexprep(curr.name,'kperm1.mat','kperm') ;
    
end

missingInd = find(diff(sort(findInds))>1)+1 ;
missingAnimal = corr_xls.CorrectedAnimalName{missingInd} ;

%% add a column 

addCol = cell(length(animal_names_rev),1) ;
addCol(findInds) = fn_pre ;

corr_xls.file_prefix = addCol ;

%% delete that row

corr_xls(missingInd,:) = [] ;
size(corr_xls,1)

%% save new table

filename = [ DD_PROC '/MamiNum_corrected_proc.txt' ] ;
writetable(corr_xls,filename,'Delimiter','\t')

