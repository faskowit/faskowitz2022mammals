clear all
close all
clc

load nets
mats = cat(3,Cij,Cij_rmvLong,Cij_rmvShort);
cn = zeros(1,3);
for i = 1:size(mats,3)
%     mats(:,:,i) = mats(:,:,i)/sum(sum(mats(:,:,i)));
%     cn(i) = f_CN_corr_aprox(mats(:,:,i));
    cn(i) = f_CN_smp(mats(:,:,i));
end