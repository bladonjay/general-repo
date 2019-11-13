function  plotRateMap(self,varargin)
%filenameps = 'C:\PlexonData\Matlab Plots\test25.ps';
filenamejpeg = 'C:\PlexonData\Matlab Plots\test32.jpeg';
X = self.x; Y = self.y;

if self.raw_pos
    self = self.FixPos;
end

if length(varargin)~=0
    cel = varargin{1} ;
else
    cel = self.cells;
    
end

xbins = min(self.b_x):10:max(self.b_x);
ybins = min(self.b_y):10:max(self.b_y);
thres = .06;
occ = histcn([X Y],xbins,ybins)*mode(diff(self.ts));
k = fspecial('gaussian',10,2);

for i = 1:size(cel,1)
    
    spks = [self.spk_x(cel(i,:))  self.spk_y(cel(i,:))];
    
    binspks = histcn(spks,xbins,ybins);
    binspks(occ==0) = nan;
    occ(occ==0) = nan;
    binspks = binspks./occ;
    binspks = nanconvn(binspks, k,'nanout',true);
    binspks(occ < thres) = nan;
    
    binspks = binspks';
    
    
    
    [im, clim] = nantowhite(binspks);
    h = figure('visible','off');
    imagesc(im,clim);
    axis off
    set(gca,'YDir','normal')
    colorbar
    title (['Cell ' num2str(self.cells(i,:))])
    %print(gcf, '-dpsc2',filenameps ,'-append')
    saveas(gcf, filenamejpeg,'jpg')
    close all
end
disp(['Saved in: ' filenamejpeg])
end