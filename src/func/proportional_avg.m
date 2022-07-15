function [pavg,w] = proportional_avg(M)          
s = sum(M,3) ; 
w = M ./ s ;
pavg = sum(M .* w, 3,'omitnan') ;
pavg(isnan(pavg)) = 0 ;