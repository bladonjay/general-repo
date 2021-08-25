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
if ~exist('F','var')
    F=nanmean(f_cond);
elseif isempty(F)
    F=nanmean(f_cond);
end

% info per pix
I=p.*(f_cond./F).*log2(f_cond./F);
% sparsity per pix
SP=(p.*f_cond.^2); % denominator
% if there are any bins with no spikes, your info for that is 0
I_nan= isnan(I);
I(I_nan)=0;
% if there are any bins with no occupancy your info for that is 0
I_zero= I==-inf;
I(I_zero)=0;
% e.g. 
Sparsity=sum(p.*f_cond)^2/sum(SP);
Skaggs_Info=nansum(I);
% JH Bladon

end



%{
% this is test code to show how the function works
probs=ones(1,50)/50; even probability of occupying each bin
cell(1,:)=normpdf(1:50,10,2)*10; % just generate some gaussian fields
cell(2,:)=normpdf(1:50,10,8)*10;
cell(3,:)=normpdf(1:50,10,20)*10;
% plot each field
figure; plot(cell'); legend('curve 1','curve 2','curve 3');
for i=1:3
[a,b]=Skaggs_basic(probs,cell(i,:),nanmean(cell(i,:)));
fprintf('for curve %d, info=%.2f, sparsity=%.2f \n',i,a,b); 
% print theirinfo and sparsity
end

%}