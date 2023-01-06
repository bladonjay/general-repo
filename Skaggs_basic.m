function [Skaggs_Info,Sparsity]= Skaggs_basic(p, f_cond, F)
% [Skaggs_Info,Sparsity]= Skaggs_basic(p, f_cond, F)

% Information: Info per second=sum over bins(p(bin).*(Ratei/mean rate).*log2(ratei/mean rate))
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

% info per spike
% I=p.*(f_cond./F).*log2(f_cond./F);
% info per second
I=p.*(f_cond).*log2(f_cond./F);

% sparsity per pix
SP=(p.*f_cond).^2; % denominator
% if there are any bins with no spikes, your info for that is 0
I_nan= isnan(I);
I(I_nan)=0;
% if there are any bins with no occupancy your info for that is 0
I_zero= I==-inf;
I(I_zero)=0;
% e.g. 
Sparsity=sum(p.*f_cond)^2/sum(SP);
Skaggs_Info=nansum(I);
% JH Bladon, originally from rob komorowski

end



%{
% this is test code to show how the function works
nbins=200;
probs=ones(1,nbins)/nbins; %even probability of occupying each bin
cell=normpdf(1:nbins,nbins*.2,nbins*.05)*10; % just generate some gaussian fields
cell(2,:)=normpdf(1:nbins,nbins*.2,nbins*.01)*10;
cell(3,:)=2-cell(1,:);

cell(4,:)=cell(1,:)+2;
cell(5,:)=5+zscore(cell(1,:));
cell(6,:)=5+-zscore(cell(1,:));

% cell(1,:) = ones(50,1);
% cell(2,:) = ones(50,1)+5;
% cell(3,:) = ones(50,1);
% cell(1,1) = 5;
% cell(2,1) = 1;
% cell(3,1) = 10;


% plot each field
figure; plot(cell'); legend('curve 1','curve 2','curve 3');
for i=1:size(cell,1)
[a,b]=Skaggs_basic(probs,cell(i,:),nanmean(cell(i,:)));
fprintf('for curve %d, info=%.8f, sparsity=%.2f \n',i,a,b); 
% print theirinfo and sparsity
end

%}