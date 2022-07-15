function [ motifM, aMat, cMat, pMat, dMat, odMat ] = comm_motif_und(cij,ca,nullModel)
% function to get community motifs, assortative, core-periphery, and
% dissassotative, as in https://www.nature.com/articles/s41467-017-02681-z
%
% INPUTS:
%       cij         original data (nans will be converted to zeros)
%       ca          |nodes| x 1 community affiliation vector
%
% OUPPUTS:
%       motifM:     matrix of motif for each off-diagonal community
%                   interaction, 1=assort, 2=core, 3=per, 4=dissort
%       aMat        edges that take part in assortative relationships
%       cMat        edges that take part in core relationships
%       pMat        edges that take part in periphery relationships
%       dMat        edges that take part in disassortative relationships
%       odMat       edges that take part in within-community relationships
%
% N.B (aMat + cMat + pMat + dMat + odMat) = cij

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    nullModel = 0 ;
end

% make sure there are no na's
cij(isnan(cij)) = 0 ;

% ensure und
cij = triu(cij,1) + triu(cij,1)' ;

% make the block mat here
[~,avgBlMat] = comm_mat_und(cij,ca) ;

k = size(avgBlMat,1) ;
n = size(cij,1) ;
motifM = zeros([k k]) ;
aMat = zeros([n n]);
cMat = zeros([n n]);
pMat = zeros([n n]);
dMat = zeros([n n]);
odMat = zeros([n n]) ;

% if we want null com motifs, randomize avgBlMat here
if nullModel == 1
    randp1 = randperm(k) ;
    randp2 = randperm(k) ;
    avgBlMat = avgBlMat(randp1,randp2) ;
end

for idx = 1:k    
    for jdx = 1:k
        
        % on diagonal
        if jdx == idx
            odMat = add_blockblock(cij,ca,jdx,idx) + odMat ;
            continue
        end
                
        % off diagonal 

        minWithinCom = min(avgBlMat(idx,idx),avgBlMat(jdx,jdx)) ;
        maxWithinCom = max(avgBlMat(idx,idx),avgBlMat(jdx,jdx)) ;
        
        subMat = [ avgBlMat(idx,idx) avgBlMat(idx,jdx) ;    % when unrolling
                   avgBlMat(jdx,idx) avgBlMat(jdx,jdx) ] ;  % 1 3 
                                                            % 2 4
        
        % position 2 comparison
        if minWithinCom > subMat(2) % assort
            motifM(jdx,idx) = 1 ;
            aMat = add_blockblock(cij,ca,jdx,idx) + aMat;
        elseif subMat(2) > maxWithinCom % dissort
            motifM(jdx,idx) = 4 ;
            dMat = add_blockblock(cij,ca,jdx,idx) + dMat ;
        else
            tmpDiff = diff([minWithinCom subMat(2) maxWithinCom]) ;
            if tmpDiff(1) < tmpDiff(2)
                % periphery because closer to smaller value
                motifM(jdx,idx) = 3 ;
                pMat = add_blockblock(cij,ca,jdx,idx) + pMat ;
            else
                % core because closer to larger value
                motifM(jdx,idx) = 2 ;
                cMat = add_blockblock(cij,ca,jdx,idx) + cMat ;
            end
        end
                             
    end  
end

function retMat = add_blockblock(dat,coms,com1,com2)
    tmp = dat(coms==com1,coms==com2) ;
    retMat = zeros(size(dat)) ;
    retMat(coms==com1,coms==com2) = tmp ;
%     retMat = mat + tmp2 ;
end

end % end func 