% fresh start
clearvars
close all

%% run config

config_file='config_mammals_1.m';
addpath(strcat(pwd,'/config'))
run(config_file);

%% load up 

filename = [ DD_PROC '/MamiNum_corrected_proc.txt' ] ;
ssheet = readtable(filename,'Delimiter','\t') ;

%%

species_names = readtable([ DD_RAW '/mami_species_names_v1.csv'],...
    'ReadVariableNames',0,'Delimiter',',') ;

db = readtable([ DD_RAW 'vertlife_taxonomies.csv' ], ...
    'Delimiter',',') ;

%%

sp_names = species_names.Var1 ;
sp_names = lower(sp_names) ;

db_names = db.scientificname ;
db_names = lower(db_names) ;

%%
findind = zeros(length(sp_names),1) ;

for x = 1:length(sp_names)
    
    disp(x)
    
    currname = sp_names{x} ;
    
    ff = find(strcmpi(currname,db_names))

    if ~isempty(ff)
        findind(x) = ff ;
    else
        findind(x) = 0 ;
    end
        
end

sp_names2 = sp_names(findind~=0) ;

%%  

sp_names2 = unique(sp_names2) ;

% https://www.mathworks.com/matlabcentral/answers/107307-function-to-capitalize-first-letter-in-each-word-in-string-but-forces-all-other-letters-to-be-lowerc
expression = '(^|\.)\s*.';
replace = '${upper($0)}';
sp_names2 = regexprep(sp_names2,expression,replace) ;

filename = [ DD_PROC '/unique_species_names_proc.csv'];
writecell(sp_names2,filename,'Delimiter',',')
