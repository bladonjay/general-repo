function [Skaggs_Info,Sparsity]= Skaggs_basic(p, f_cond, F)
% [Skaggs_Info,Sparsity]= Skaggs_basic(p, f_cond, F)

% Information: Info=sum over bins(p(bin).*(Ratei/mean rate).*log2(ratei/mean rate))
% Sparsity: Sparsity= sum over i bins(Pi*Ratei^2) / mean rate^2
% p is the occupancy of that pixel (prob of occupancy)
% F_cond= rate at that pixel
% F = average firing rate

% any bins that never got visited?
novisits=p==0;
p(novisits)=[]; f_cond(novisits)=[];

% info per pix
I=p.*(f_cond./F).*log2(f_cond./F);
% sparsity per pix
SP=(p.*f_cond.^2);
% if there are any bins with no spikes, your info for that is 0
I_nan= isnan(I);
I(I_nan)=0;
% if there are any bins with no occupancy your info for that is 0
I_zero= I==-inf;
I(I_zero)=0;
Sparsity=sum(SP)/F^2;
Skaggs_Info=nansum(I);
% JH Bladon

end