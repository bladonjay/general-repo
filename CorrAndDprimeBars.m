function [hfig, allbars,legend] = CorrAndDprimeBars(cor,whichplot)
%Builds Correlation and Dprime bars item*pos, item, pos and context for the correlations
%sam made. cor is your cor from sams files, whichplot is either 1,2 or 3. 1
%is just correlation, 2 is just dprime and 3 is both

if ~exist('whichplot','var')
    whichplot=3;
elseif isempty(whichplot)
    whichplot=3;
end
    



CPI=cell2mat(arrayfun(@(a) a.cor.CPIm(:),cor,'uni',0)');
CPxI=cell2mat(arrayfun(@(a) a.cor.CPxIm(:),cor,'uni',0)');
CP=cell2mat(arrayfun(@(a) a.cor.CPm(:),cor,'uni',0)');
CxP=cell2mat(arrayfun(@(a) a.cor.CxPm(:),cor,'uni',0)');
CxPI=cell2mat(arrayfun(@(a) a.cor.CxPIm(:),cor,'uni',0)');
CxPxI=cell2mat(arrayfun(@(a) a.cor.CxPxIm(:),cor,'uni',0)');
xC=cell2mat(arrayfun(@(a) a.cor.xCm(:),cor,'uni',0)');
xCxPI=cell2mat(arrayfun(@(a) nanmean(a.cor.xCxPIm(:)),cor,'uni',0)');
xCxPxI=cell2mat(arrayfun(@(a) nanmean(a.cor.xCxPxIm(:)),cor,'uni',0)');

match=nanmean([CPI CPxI CxPI CxPxI  xCxPI xCxPxI],1);
matchSEM=SEM([CPI CPxI CxPI CxPxI   xCxPI xCxPxI],1);
allbars.rawcor=[CPI CPxI CxPI CxPxI  xCxPI xCxPxI];
allbars.cor=match;
allbars.cor(2,:)=matchSEM;
% always make a figure
if any(whichplot==[1 2 3])
    hfig=figure;
    
else
    hfig=0;
end
if whichplot==1 || whichplot==3
    if whichplot==3
        subplot(1,2,1);
    end
    
    bar(match,.5,'b'); hold on
    hh = errorbar(match,matchSEM,'linestyle','none','color','k');
    set(gca,'XTick',1:6,'XTickLabel',{'CPI','CPxI','CxPI','CxPxI','xCxPI','xCxPxI'},'fontsize',14)
    %xticklabel_rotate([],45);
    ylabel('Mean Corr. Coef.')
    
end
%%

itemXposR=cell2mat(arrayfun(@(a) dprime(a.cor.CPI(:),a.cor.CPxI(:)),cor,'uni',0)');
posR=cell2mat(arrayfun(@(a) dprime(a.cor.CP(:),a.cor.CxP(:)),cor,'uni',0)');
itemR=cell2mat(arrayfun(@(a) dprime(a.cor.CxPI(:),a.cor.CxPxI(:)),cor,'uni',0)');
contextR=cell2mat(arrayfun(@(a) dprime(a.cor.CxP(:),a.cor.xC(:)),cor,'uni',0)');

valenceR=cell2mat(arrayfun(@(a) dprime([a.cor.CI(:); a.cor.xCxPxI(:)],[a.cor.CxI(:);a.cor.xCI]),cor,'uni',0)');
itemXcontextR=cell2mat(arrayfun(@(a) dprime(a.cor.CxPI(:),a.cor.xCxPxI(:)),cor,'uni',0)');

ratio=nanmean([contextR posR itemR itemXposR valenceR itemXcontextR]);
ratioSEM=SEM([contextR posR itemR itemXposR valenceR itemXcontextR]);
allbars.rawdprime=[contextR posR itemR itemXposR valenceR itemXcontextR];
allbars.dprime=ratio;
allbars.dprime(2,:)=ratioSEM;

if whichplot==2 || whichplot==3
    if whichplot==3
        subplot(1,2,2);
    end
    bar(ratio'); hold on;
    
    barcenters(:,1)=[1:6];
    hh = errorbar(barcenters,ratio',ratioSEM','linestyle','none','color','k');
    set(gca,'XTickLabel',{'Context','Pos','item','item/pos','Valence','Item*Context'})
    ylabel('Dprime Magnitude');
end

legend={'CPI' 'CPxI' 'CxPI' 'CxPxI'  'xCxPI' 'xCxPxI';...
    'Context' 'Pos' 'item' 'item/pos' 'Valence' 'Item*Context'};



end

