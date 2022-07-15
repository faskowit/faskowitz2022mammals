function [WMD,mf] = f_interpolate_WMD(L,ED)

N = size(L,1);

% correct for fibers that are shorter than the ED
ff = find((L>0)&((L-ED)<0));
minED = min(min(ED + (1000*eye(N))));
L2remove = find((ED(ff) - L(ff)) > minED);
L(ff(L2remove)) = 0;

ffmask = find(triu(ones(N),1));
[~,orderEDind] = sort(ED(ffmask),'ascend');
rowk = zeros(1,length(orderEDind)); colk = zeros(1,length(orderEDind));
for k=1:length(orderEDind)
    M = zeros(N);
    M(ffmask(orderEDind(k))) = 1;
    [rowk(k) colk(k)] = find(M);
end;

ff = find(L>0);

% linear interpolation
%lin_pol = polyfit(ED(ff), L(ff), 1);
%lin_vals = polyval(lin_pol, ED(ff));

%  second ord. interpolation
sec_pol = polyfit(ED(ff), L(ff), 2);
sec_vals = polyval(sec_pol, ED(ff));

% figure; plot(ED(ff), L(ff), '.');
% hold on
% plot(ED(ff), sec_vals, 'g.');
% hold off;

L0 = L;

thr = 0; Lthr = zeros(N); mf = zeros(N);
rad_max = .2;

% interpolate in random order
lk = length(rowk);
rp = randperm(lk);
rowk = rowk(rp); colk = colk(rp);

for k=1:length(rowk)
    i = rowk(k); j = colk(k);
%    thr = 0; maxthr = ceil(rad_max*ED(i,j));
    thr = 0; maxthr = max(10,ceil(rad_max*ED(i,j)));
    flag_maxthr = 0;
    while (L(i,j)==0) && ~flag_maxthr
        thr = thr+1;
        if thr == maxthr
            L(i,j) = polyval(sec_pol, ED(i,j));
            flag_maxthr = 1;
        else
            neiI = find(ED(i,:)<thr);
            neiJ = find(ED(j,:)<thr);
            if (isempty(neiI))&&(isempty(neiJ))
                L(i,j) = 0;
            else
                dj_neiI = 0; di_neiJ = 0; faci_neiJ = 0; facj_neiI = 0;
                if any(L(j,neiI))
                    dj_neiI = mean(nonzeros(L(j,neiI)));
                    facj_neiI = L(j,neiI)./ED(j,neiI);
                end
                if any(L(i,neiJ))
                    di_neiJ = mean(nonzeros(L(i,neiJ)));
                    faci_neiJ = L(i,neiJ)./ED(i,neiJ);
                end
                if any([di_neiJ, dj_neiI])
                    %dIJ = mean(nonzeros([di_neiJ, dj_neiI]));
                    %L(i,j) = max(polyval(sec_pol, ED(i,j)), dIJ);
                    meanfac = mean(nonzeros([facj_neiI, faci_neiJ]));
                    if (meanfac>4) meanfac = 4; end;
                    mf(i,j) = meanfac;
                    L(i,j) = meanfac.*ED(i,j);
                    %L(i,j) = dIJ;
                else
                    L(i, j) = 0;
                end
            end
        end
    end
    Lthr(i,j) = thr;
end;

Lsec = triu(L,1);
Lsec = Lsec+Lsec';
%figure; plot(ED, Lsec, '.');

WMD = Lsec;

