function DI=dunns(distM,ca)   
%%%Dunn's index for clustering compactness and separation measurement
% dunns(clusters_number,distM,ind)
% clusters_number = Number of clusters 
% distM = Dissimilarity matrix
% ind   = Indexes for each data point aka cluster to which each data point
% belongs
uu=unique(ca) ;
nn=length(uu);

denominator=[];
for idx=1:nn
    indi=find(ca==uu(idx));
    indj=find(ca~=uu(idx));
    x=indi;
    y=indj;
    temp=distM(x,y);
    denominator=[denominator;temp(:)];
end
num=min(min(denominator)); 
neg_obs=zeros(size(distM,1),size(distM,2));
for ix=1:nn
    indxs=find(ca==ix);
    neg_obs(indxs,indxs)=1;
end
dem=neg_obs.*distM;
dem=max(max(dem));
DI=num/dem;
end

%% 

for idx=1:nn
    indi=find(ca==uu(idx));
    indj=find(ca~=uu(idx));
    x=indi;
    y=indj;
    temp=distM(x,y);
    denominator=[denominator;temp(:)];
end



