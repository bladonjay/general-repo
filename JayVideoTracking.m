function ProcessedTracking = JayVideoTracking(vidfilepath,coorddata,arenamask)

% ProcessedTracking = JayVideoTracking(vidfilepath,coorddata,arenamask)
% vidfilepath is an avi file
% coorddata is nframes by at least 5 col matrix. col 1 is timestamp of each
% frame, this is crucial.  The next 4 columns are the x and y of the old
% tracking

% most important part here is the cutoffs;
% line 109 is a 90% of max pixel, this is hardcoded but you may want to
% change it
% line 127,8,9 are where the color cutoffs are.  You want a mean luminance
% thats high, and you want a colorness thats high as well.  You'll want to
% change these based on the video. If the video is saturated, you'll have
% trouble digging out your colors because it will all be white.
% JHB 12/21/18

% use this to get a handle on the average luminance of the LED's
tryoldtracking=0;


vidfile=VideoReader(vidfilepath);
frame1=read(vidfile,1);
% if any dimension of our arenamask isnt the same size as our video file
if any(size(arenamask)~=size(frame1(:,:,1)))
    [newframe3]=brightengui(uint8(frame1),'bottom');
    [arenamask,ix,iy]=roipoly(newframe3);
    %SuperRat(k).delaydoor=arenamask;
    close(gcf);
end


% this fx takes way too long
%vidinfo=aviinfo(SuperRat(i).vidfile);
vidfile=VideoReader(vidfilepath);

nframes=round(vidfile.Duration*vidfile.FrameRate);
% read the first 2 minutes seconds and cat those into an average luminance;
allframe=[];
% 1000 frames over about 2 minutes
for i=1:4:4000
    % grab the frame
    [frame1]=double(readFrame(vidfile))/1000;
    % each frame comes in as int8, but it needs to be converted to
    % uint8. thus, it needs to be brightened substantially
    
    if i==1, allframe=frame1;else, allframe=[allframe+frame1]; end
end
avgframe=allframe;
% this gets us a good background, we could keep updating our background
% as we go, just in case it gets super bright... we could also maybe
% put in a stop if that happens

% now I think this should be the training dataset here, so lets check
% the real tracking

% grab arena mask
vidfile=VideoReader(vidfilepath);
frame1=read(vidfile,1);
%coorddata=SuperRat(k).coords;
coordname=fieldnames(coorddata);
coords=coorddata.(coordname{1});
% doctor coords so theyre only 5 columns
coords(:,sum(coords)==0)=[];
coords(coords==0)=nan;
% brighten image so we can see it
%[newframe3]=brightengui(uint8(frame1),'bottom');
%arenamask=SuperRat(k).arenamask;
% now that we have the mask we can run the session;

%figure('position',[156 85 690 880]);
goodcoord=[];
h=waitbar(0,'testing');
maxframe=min([nframes size(coords,1)]) ;
mycoords=[];
for i=1:2:maxframe
    % every say 2 minutes, or 900 frames (or if were doing every frame every 3200)
    % recalculate the background image
    if mod(i,3600)==0
        allframe=[];
        % weight it forward so that we catch the average frame
        % prospectively
        for bg=-1500:6:4500
            [frame1]=double(read(vidfile,i+bg))/1000;
            if isempty(allframe)
                allframe=frame1;
            else
                allframe=[allframe+frame1];
            end
        end
    end
    frame=read(vidfile,i);
    %subplot(2,1,1);
    %image(flipud(int8(frame)));
    %hold on;
    %set(gca,'YDir','normal')
    if tryoldtracking
        if ~isnan(coords(i,2)) && ~isnan(coords(i,3)) % both colors
            %plot(nanmean(coords(i,[2 4]))*.6246,nanmean(coords(i,[3 5]))*.6246,...
            %    'w+','MarkerSize',12);
        end
        if ~isnan(coords(i,2)) %
            % plot(coords(i,2)*.6246, coords(i,3)*.6246,...
            %     'r+','MarkerSize',12);
            mycolors=frame(round(coords(i,3)*.6246),round(coords(i,2)*.6246),:);
            coords(i,6)=nanmean(mycolors); coords(i,7)=mycolors(1)-nanmean(mycolors(2:3));
        end
        if ~isnan(coords(i,4))
            %  plot(coords(i,4)*.6246, coords(i,5)*.6246,...
            %     'g+','MarkerSize',12);
            mycolors=frame(round(coords(i,5)*.6246),round(coords(i,4)*.6246),:);
            coords(i,8)=nanmean(mycolors); coords(i,9)=mycolors(2)-nanmean(mycolors(1:3));
        end
    end
    %hold off;
    %subplot(2,1,2);
    %image(int8(frame));
    %hold on;
    subtractframe=round(nanmean(abs(double(frame)-avgframe),3),2);
    subtractframe(~arenamask)=0;
    [a,b,c]=unique(linearize(subtractframe(arenamask)));
    cutoff=prctile(a,90);
    highframe=subtractframe; highframe(highframe<cutoff)=nan;
    temp=regionprops(subtractframe>cutoff,'Centroid','Area','PixelIdxList');
    for b=1:length(temp)
        temp(b).meanheat=nanmean(subtractframe(temp(b).PixelIdxList));
        % rgb here
        for j=1:3
            % grab this color
            colorim=frame(:,:,j);
            % get mean for that color
            temp(b).meancolor(j)=nanmean(colorim(temp(b).PixelIdxList));
        end
        temp(b).redness=temp(b).meancolor(1)-nanmean(temp(b).meancolor([2 3]));
        temp(b).greenness=temp(b).meancolor(2)-nanmean(temp(b).meancolor([1 3]));
    end
    % now red ought to be whichever beats thresholds
    % brighness thresh will be above
    structcell=struct2cell(temp); cellmat=cell2mat(structcell([4 6 7],:));
    bigdots=cellmat(1,:)>30;
    reddots=cellmat(2,:)>20;
    greendots=cellmat(3,:)>20;
    % first get our thresholded
    if any(bigdots & (reddots | greendots))
        % first red
        if any(bigdots & reddots)
            redidx=find(bigdots & reddots,1,'first');
            mycoords(i,2:3)=temp(redidx).Centroid;
            % plot(mycoords(i,2), mycoords(i,3),...
            %'r+','MarkerSize',12);
            mycoords(i,6)=nanmean(temp(redidx).meancolor);
            mycoords(i,7)=temp(redidx).meancolor(1)-nanmean(temp(redidx).meancolor(2:3));
        end
        % then green
        if any(bigdots & greendots)
            greenidx=find(bigdots & greendots,1,'first');
            mycoords(i,4:5)=temp(greenidx).Centroid;
            % plot(mycoords(i,4), mycoords(i,5),...
            %'g+','MarkerSize',12);
            mycoords(i,8)=nanmean(temp(greenidx).meancolor);
            mycoords(i,9)=temp(greenidx).meancolor(2)-nanmean(temp(greenidx).meancolor(1:3));
        end
    end
    %hold off;
    %pause(.03);
    waitbar(i/maxframe,h,sprintf('still testing rat %s', 'a rat...'));
end
waitbar(1,h,'done testing');
pause(.1); close(h);
ProcessedTracking.coords=coords;
ProcessedTracking.mycoords=mycoords;

%% this fills in the remaining timestamps and gets your average rat position
%
%
%
%
%
%
summarycoordds=[];
mycoords=ProcessedTracking.mycoords;
coordts=ProcessedTracking.coords(:,1);
mycoords(:,1)=coordts(1:length(mycoords));
mycoords(mycoords==0)=nan;
% interpolate everything thats flanked by two real numbers
% gonna do a linear interp
for i=2:length(mycoords)-1
    for j=2:9
        if isnan(mycoords(i,j))
            mycoords(i,j)=mean(mycoords([i-1 i+1],j));
        end
    end
end
% now pick green, but if its empty go for red. red is bad becuase there
% is an indicator light on the maze that is bad for study period
summarycoords=mycoords(:,[1 4 5]);
for i=1:length(summarycoords)
    if isnan(summarycoords(i,2))
        summarycoords(i,[2 3])=mycoords(i,[2 3]);
    end
end
ProcessedTracking.summarycoords=summarycoords;
fprintf('Session done \n');






end

