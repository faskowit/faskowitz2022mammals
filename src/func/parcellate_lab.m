function [outVec, outVec2] = parcellate_lab(dat,labVec,labs,funcopt)

if length(dat) ~= length(labVec)
   error('data and label need to be same length') 
end

if nargin < 3 
   labs = unique(labVec) ; 
end

if nargin < 4 
   funcopt = 'mean' ; 
end

% ensure column
if ~iscolumn(dat) || ~iscolumn(labVec) || ~iscolumn(labs) 
   error('inputs should be column arrays')  
end

%% run it

switch funcopt
    case {'','mean'}
        reducefunc = @mean ;
    case {'sum'}
        reducefunc = @sum ;
    case {'count'}
        reducefunc = @length ;
    case {'std'}
        reducefunc = @std ;
    case {'nanmean'}
        reducefunc = @nanmean ;
    case {'nansum'}
        reducefunc = @nansum ;
    otherwise
        error([ 'not a viable option for func ' funcopt ])
end
        
% outVec = zeros(length(labs),1) ;
% for idx = 1:length(labs)
%     currLab = labs(idx) ;
%     outVec(idx) = mean(dat(labVec==currLab));
% end

outVec = arrayfun(@(x) reducefunc(dat(labVec==x)), labs) ;

if nargout > 1
%     % vector in original length
%     outVec2 = zeros(size(dat)) ;
%     for idx = 1:length(labs)
%         currLab = labs(idx) ;
%         outVec2(labVec==currLab) = outVec(idx);
%     end
    dumArray = labVec == labs' ;
    [dumInds,~] = find(dumArray') ;
    outVec2 = outVec(dumInds) ; 
end


