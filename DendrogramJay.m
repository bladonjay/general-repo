function [hfig,hfig2,H1, H2] = DendrogramJay(cor,notevenodd)
% function DendrogramJay(corr, not even or odd)
%
%Makes a Dendrogram figure, with the handle hfig
% uses the average value for each cluster (trial type), and calculates
% the correlation between each.  The cor field has to have a 16 row matrix
%with at least 5 columns in the cor.cor(i).spk1 matrix. There should be 8
%different trial types each cut in half.

% outputs
%   hfig= handle for each figure. if you dont put a 1 or a 2 for evenodd,
%   hfig2 is just 0
%   H1,H2=handles for the dendrogram bars to change color or weight

% notevenodd is the single type breakdown.
if ~exist('notevenodd','var')
    notevenodd=0;
elseif isempty(notevenodd)
    notevenodd=0;
end

rate=cell2mat(arrayfun(@(a) a.cor.spk1',cor,'uni',0)');

rate = rate(~any(isnan(rate),2),:)';
Z1 = linkage(rate,'average','correlation');
hfig=figure;
[H1,T,outperm]=dendrogram(Z1,16);

% names1 = {'A1+','A1+','B1-','B1-','A2+','A2+','B2-','B2-'...
%     'A3-','A3-','B3+','B3+','A4-','A4-','B4+','B4+'};
names1 = {'X1+','X1+','Y1-','Y1-','X2+','X2+','Y2-','Y2-'...
    'X3-','X3-','Y3+','Y3+','X4-','X4-','Y4+','Y4+'};
standardMat=[];
set(gca,'XTickLabel',names1(outperm),'fontsize',13)

set(H1([1:8]),'linewidth',2) % set(H1([1:7,9]),'linewidth',2)
set(H1([9:12]),'linewidth',4) % set(H1([8,10:12]),'color','g','linewidth',4)
set(H1([14 13]),'linewidth',5) % set(H1([14 13]),'color','r')
set(H1([15]),'linewidth',6) % set(H1([15]),'color','k')
set(gca,'YTick',.1:.1:1.3, 'YTickLabel',1-(.1:.1:1.3))


ylabel('Mean Corr. Coef')

if notevenodd==1
    names1 = {'A1+','B1-','A2+','B2-',...
        'A3-','B3+','A4-','B4+'};
    hfig2=figure;
    standardMat=[];
    rate=cell2mat(arrayfun(@(a) a.cor.spk',cor,'uni',0)');
    rate = rate(~any(isnan(rate),2),:)';
    Z1 = linkage(rate,'average','correlation');
    [H2,T,outperm]=dendrogram(Z1,8);
    set(gca,'XTickLabel',names1(outperm),'fontsize',13)
    %
    set(H2([1:4]),'linewidth',2) %set(H2([1:4]),'color','g')
    set(H2([5 6]),'linewidth',4) %set(H2([5 6]),'color','r')
    set(H2([7]),'linewidth',5) %set(H2([7]),'color','k')
    set(gca,'YTick',.1:.1:1.3, 'YTickLabel',1-(.1:.1:1.3))
    
    
    ylabel('Mean Corr. Coef')
elseif notevenodd==2
    
    % this calculates the dendrogram based on the averages by half, and
    % should sort out exactly like the even odd one wth 16 categories
    rate=cell2mat(arrayfun(@(a) a.cor.spk1',cor,'uni',0)');
    rate = rate(~any(isnan(rate),2),:)';
    Z1 = linkage(rate,'average','correlation');
    hfig2=figure;
    % T is where each single name sorted out to
    
    [H2,T,outperm]=dendrogram(Z1,8);
    names1 = {'A1+','A1+','B1-','B1-','A2+','A2+','B2-','B2-'...
        'A3-','A3-','B3+','B3+','A4-','A4-','B4+','B4+'};
    % lets check that the nodes paired out to
    
    for i=1:length(outperm)
        
        temp=find(T==outperm(i),2);
        if length(temp)>1
            if strcmpi(names1{temp(1)},names1{temp(2)})
                ordered{i}=names1{temp(1)};
            else
                ordered{i}=[names1{temp(1)} ' / ' names1{temp(2)}];
            end
        else
            ordered{i}=[names1{temp(1)} ' / ?'];
        end
    end

    set(gca,'XTickLabel',ordered,'fontsize',13)
    
    set(H2([1:4]),'linewidth',2)
    set(H2([5:6]),'linewidth',4)
    set(H2([7]),'linewidth',5)
    
    set(gca,'YTick',.1:.1:1.3, 'YTickLabel',1-(.1:.1:1.3))
    
    
    ylabel('Mean Corr. Coef')
    
    

else
    hfig2=0;
end


end

