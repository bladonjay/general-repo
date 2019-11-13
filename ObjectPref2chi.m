function [chivals,pobj,ratios] = ObjectPref2chi(session,unitdata,obj1coords,obj2coords)
% function [chivals]=ObjectPref2chi(session,unitdata,sampledist)
% this function runs a chi square test on the occupancy and spike rates of
% a rat near my two objects.  Its mostly still hardcoded but i'm working on
% it.
% Session needs fields session.edit_coords, and unitdata has to have field
% unitdata.units(:).ts
sampledist=40; chivals=[]; suppress=0;

% first coords
coords=session.edit_coords;

coords(:,2)=smooth(coords(:,2)); coords(:,3)=smooth(coords(:,3));
displacement=sqrt(diff(coords(:,2)).^2+diff(coords(:,3)).^2);
% 30 frames/sec and about .6 cm/pixel to go from pix/frame to
% cm/sec.  I'm not gonna threshold because sampling occurs at a standstill
velocity=[1; (displacement.*[diff(coords(:,1))])*20]; % 30 frames/sec
if ~exist('obj1coords','var')
    obj1coords=[501,599];
end
if ~exist('obj2coords','var')
    obj2coords=[301,236];
end

clear obj1time obj2time
% grab the time he spends close to the objects
obj1time(:,1)=coords(:,2)-obj1coords(1);
obj1time(:,2)=coords(:,3)-obj1coords(2);
obj1tdist=sqrt(obj1time(:,1).^2+obj1time(:,2).^2);
obj1sampleframes=obj1tdist<sampledist;  % 3 pixels, or 15 cm
obj1sampletime=sum(obj1sampleframes*.033);

obj2time(:,1)=coords(:,2)-obj2coords(:,1);
obj2time(:,2)=coords(:,3)-obj2coords(:,2);
obj2tdist=sqrt(obj2time(:,1).^2+obj2time(:,2).^2);
obj2sampleframes=obj2tdist<sampledist;
obj2sampletime=sum(obj2sampleframes*.033);

awayfromobjects=~(obj1sampleframes & obj2sampleframes);
% coords2=coords(velocity>0.5,:);
for i=1:length(unitdata.units)
    cellspikes=unitdata.units(i).ts;
    okspikes=cellspikes(cellspikes>min(coords(:,1)) & cellspikes<max(coords(:,1)));
    if length(okspikes)< 24000 && length(okspikes)>80 % average above 40 hz?
        % snap them to pixels and timestamps
        
        sesspikes = interp1(coords(:,1),coords(:,1:end),okspikes,'nearest');
        if suppress
            figure;
            plot(coords(:,2), coords(:,3), 'k*');
            hold on;
            %rat path for entire session, regardless of velocity
            
            plot(sesspikes(:,2), sesspikes(:,3), 'r.');
            hold off
        end
        % get spikes for first object
        clear obj1spikes obj2spikes
        % get the distance of each spike from the object
        obj1spikes(:,1)=sesspikes(:,2)-obj1coords(:,1);
        obj1spikes(:,2)=sesspikes(:,3)-obj1coords(:,2);
        % now take all the timestamps he was close to the object
        obj1dist=sqrt(obj1spikes(:,1).^2+obj1spikes(:,2).^2);
        % now only take the spikes that are close to that object
        
        obj1spiketime=obj1dist<sampledist;
        % get spikes for second object
        obj2spikes(:,1)=sesspikes(:,2)-obj2coords(:,1);
        obj2spikes(:,2)=sesspikes(:,3)-obj2coords(:,2);
        obj2dist=sqrt(obj2spikes(:,1).^2+obj2spikes(:,2).^2);
        
        % now how long was he close to the object?
        obj2spiketime=obj2dist<sampledist;
        % now compare all: top row is time spent at each object, bottom row
        % is the number of spikes that occurred at each object
        chimat=([round(obj1sampletime) round(obj2sampletime);
            sum(obj1spiketime) sum(obj2spiketime)]);
        
        [pval,chistat,df]=chi2indep(chimat);
        %chivals(i)=chi2pdf(chistat,df);
        chivals(i)=chistat;
        ratios(i,:)=(chimat(1,:)./chimat(2,:))';
        %pobj(i)=diff(ratios)/sum(ratios);
        pobj(i)=pval;
    else
        pobj(i)=nan;
        chivals(i)=1;
        ratios(i,1:2)=nan;
    end
end






end

