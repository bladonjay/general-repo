function [Info,Sparsity,InfoP,SparsityP] = Bootstrap_Skaggs(p,curve,nboots)
% DONT USE, THIS WILL NOT WORK AS OF RIGHT NOW %
% [Info,Sparsity,InfoP,SparsityP] = Bootstrap_Skaggs(p,curve,mymean,nboots)

% bootstraps skaggs information for a spiking neuron dataset.  bins are
% independent, and info is returned as bits, which you have to multiply by
% p(spike)
verbose=1;
mymean=nanmean(curve);
[Info,Sparsity]=Skaggs_basic(p,curve,mymean);

% mixup for the boots
bootinfo=nan(1,nboots); bootsparsity=nan(1,nboots);
parfor i=1:nboots
    [bootinfo(i),bootsparsity(i)]=Skaggs_basic(p,curve(randperm(length(curve))),mymean);
end
% how many boots did we beat
InfoP=nanmean(Info>bootinfo); SparsityP=nanmean(Sparsity<bootsparsity);

% plot if we ask it to
if verbose
    figure; h=histogram(bootinfo); hold on;
    line([Info Info],[0 max(h.BinCounts)],'r');
    
    figure; h=histogram(bootsparsity); hold on;
    line([Sparsity Sparsity],[0 max(h.BinCounts)],'r');
end

end

