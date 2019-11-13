function fixPosition(unitdata,Cinedata)
import CMBHOME.Utils.*
X1=cellfun(@str2num,Cinedata.markers{1}.values{1}.strings);
Y1=cellfun(@str2num,Cinedata.markers{1}.values{2}.strings);
ts1=Cinedata.markers{2}.timestamps;


X2=cellfun(@str2num,Cinedata.markers{2}.values{1}.strings);
Y2=cellfun(@str2num,Cinedata.markers{2}.values{2}.strings);
ts2=Cinedata.markers{2}.timestamps;

time=[min([ts1 ;ts2]):mode(ts1):max([ts1 ;ts2])];
time1=time;
time2=time;

X1=avghist(ts1,X1,time)';
Y1=avghist(ts1,Y1,time)';
X2=avghist(ts2,X2,time)';
Y2=avghist(ts2,Y2,time)';

X1(isnan(X1))=X2(isnan(X1));
X2(isnan(X1))=X1(isnan(X1));
Y1(isnan(X1))=Y2(isnan(X1));
Y2(isnan(X1))=Y1(isnan(X1));


jitter_threshold = 5/.25;

max_allowed_flips = 4; % samples


flips = findOnsetsAndOffsets(isnan(X1));

flips(:,2) = flips(:,2)+1;

flips = cat(1, 1, flips(:), find([0; sqrt(diff(X1).^2 + diff(Y1).^2)]>jitter_threshold), length(X1));

flips = sort(unique(flips)); % indeces of NaNs or jumps

flips = [flips(1:end-1), flips(2:end)];  % epochs formation

flips(flips(:,2)-flips(:,1)>max_allowed_flips,:) = [];

flips(:,2) = flips(:,2)-1; % adjust for diff shift

flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps

flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices

X1([flips{:}]) = nan; % remove samples in ts and X1
Y1([flips{:}]) = nan; % remove samples in ts and X1
time1([flips{:}]) = [];

X1 = interp1(time1, X1, time);
Y1 = interp1(time1, Y1, time);

X1 = ndnanfilter(X1, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
Y1 = ndnanfilter(Y1, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);

flips = findOnsetsAndOffsets(isnan(X2));

flips(:,2) = flips(:,2)+1;

flips = cat(1, 1, flips(:), find([0; sqrt(diff(X2).^2 + diff(Y2).^2)]>jitter_threshold), length(X2));

flips = sort(unique(flips)); % indeces of NaNs or jumps

flips = [flips(1:end-1), flips(2:end)];  % epochs formation

flips(flips(:,2)-flips(:,1)>max_allowed_flips,:) = [];

flips(:,2) = flips(:,2)-1; % adjust for diff shift

flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps

flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices

X2([flips{:}]) = []; % remove samples in ts and X1
Y2([flips{:}]) = []; % remove samples in ts and X1
time2([flips{:}]) = [];

X2 = interp1(time2, X2, time);
Y2 = interp1(time2, Y2, time);

X2 = ndnanfilter(X2, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
Y2 = ndnanfilter(Y2, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);


for i=1:length(initdata.unit)
    spkts=unitdata.units(i).ts;
    
    
end