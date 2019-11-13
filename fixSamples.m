function samples = fixSamples(samples,door,begin,context)
%for samples
% First column is timestamp
% Second column is correct (1) / incorrect (0)
% Third column is left (1) / right (2)
% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
% Fifth column is west (context 1) or east (context 2)
% Sixth column is position (1-4)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)
% Eighth column is Odor (A = 1, B = 2, C = 3, D = 4)
% Ninth column is duration of the sample.
% Tenth column is the day number
% Eleventh column is actual trial #
% Twelth column is the # of samples on that pot
% Thirteenth column is rat being correct (dug on cor. no dig on incorrect)
% Fourteenth column is last sample of that pot for that trial)
% Fifteenth column is whether to exclude the trial

[~,longest] = max([length(door) length(begin) length(context)]);

switch longest
    case 1
        ts=door(:);
    case 2
        ts=begin(:);
    case 3
        ts=sort(context(:));
end

[~,bin]=histc(samples(:,1),[ts;inf]);
allbins=unique(bin);
for i=1:length(allbins)
    
    tr=samples(bin==allbins(i),:);
    ind=find(bin==allbins(i));
    [o,a,b]=unique(tr(:,2:7),'stable','rows');
    if length(unique(b))>2
        disp(['Too many different conditions for the same trial see' num2str(round(10*tr(1,1))/10)])
        keyboard
    else
        
        temp=b;
        b(temp==1)=cumsum(temp(temp==1)==1);
        b(temp==2)=cumsum(temp(temp==2)==2);
        
        lastTr=zeros(length(b),1);
        lastTr(find(temp==1,1,'last'))=1;
        lastTr(find(temp==2,1,'last'))=1;
        
        
        if  tr(end,1)-tr(1,1)>30
            warning=ones(length(b),1);
        else
            warning=zeros(length(b),1);
        end
        
    end
    
    samples(ind,11)=i;
    samples(ind,12)=b;
    samples(ind( tr(:,2)'==1 & ind(1):ind(end)<ind(end)),13) = 0; %correct and left
    samples(ind( tr(:,2)'==1 & ind(1):ind(end)==ind(end)),13) = 1; %correct and dug
    samples(ind( tr(:,2)'==0 & ind(1):ind(end)==ind(end)),13) = 0; %incorrect and dug
    samples(ind( tr(:,2)'==0 & ind(1):ind(end)<ind(end)),13) = 1; %incorrect and didn't dig
    samples(ind,14)=lastTr;
    samples(ind,15)=warning;
    
    
end


    
    

end