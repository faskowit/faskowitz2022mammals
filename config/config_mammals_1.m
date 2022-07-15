%% PROJECT CONFIG

PROJ_DIR = pwd ; % put base path to project her
cd(PROJ_DIR)

%% add to the path

projPathDirs = { 
    'src' 
    'data'
    'bin'
} ;

for idx=1:length(projPathDirs)
    addpath(genpath(strcat(PROJ_DIR,'/',projPathDirs{idx})))
end

clear projPathDirs

%% setup global vars

OUTSTR = 'mamv3' ;
LOADSTR = 'run1' ;

%% make output directory vars

DATADIR = strcat(PROJ_DIR , '/data/') ; 
DD_INTERM = strcat(DATADIR, '/interim/' ) ;
DD_PROC = strcat(DATADIR, '/processed/' ) ;
DD_RAW = strcat(DATADIR, '/raw/' ) ;

%% stuff

THRDENS = [ 0 0.05 0.1 0.15 ] ;
NNODES = 200 ;
NRESAMP = 100 ; 


