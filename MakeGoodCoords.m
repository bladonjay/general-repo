function [coorddata]=MakeGoodCoords(avifilepath,sessionname)
% annotations to come, this is hot off the script
if ~exist('avifilepath','var')
    [avifile,avipath]=uigetfile('.avi','What avi file would you like','/C');
    avi_filepath=[avipath avifile];
    vidfile=VideoReader(avi_filepath);
elseif isempty(avifilepath)
    [avifile,avipath]=uigetfile('.avi','What avi file would you like','/C');
    avi_filepath=[avipath avifile];
    vidfile=VideoReader(avi_filepath);
end


[dvtfile]=uigetfile('.DVT','now choose dvt file',avipath);
dvt_filepath=[avipath dvtfile];
pos_data=importdata(dvt_filepath);

if ~exist('sessionname','var')
    sessionname='Aladdin Session 3 11-3-16';
elseif isempty(sessionname)
    sessionname='Aladdin Session 3 11-3-16';
end

framenum=100;
fprintf('outline the maze, right click to add mask \n');
loadframe=read(vidfile,framenum);
arenamask=roipoly(loadframe);

meanimage=[];
ticker=waitbar(0,'gathering background');
for i=1:5:1000
    loadframe=read(vidfile,i);
    typechg=double(loadframe)/200;
    if i==1
        meanimage=typechg;
    else
        meanimage=meanimage+typechg;
    end
    waitbar(i/1000)
end
close(ticker); clear ticker

% do every third frame, and only take pixels within 30 mm away


circleinds=[sin(-pi:.2:pi)' cos(-pi:.2:pi)']*30;

% loading large video file chunks makes this a million times faster
skipby=3;
totframes=round(vidfile.duration*vidfile.framerate);
% now all is with regard to the frames we have
videochunks=5;
% list of each frame we'll pull
framelist=1: skipby: totframes;
% and which chunk these frames will be pulled
allchunks=ceil((1:length(framelist))./(length(framelist)/videochunks));

for chunk=2:videochunks
    % first pull the frames you want
    firstidx=find(allchunks==chunk,1,'first');
    firstframe=framelist(firstidx);
    lastidx=find(allchunks==chunk,1,'last');
    lastframe=framelist(lastidx);
    
    % have to load all frames even ones inbetween
    frameload=read(vidfile,[firstframe lastframe]);
    
    % now go from first frame of this chunk to last and read
    for k=framelist(allchunks==chunk)
        
        
        loadframe=frameload(:,:,:,((k-firstframe)/3)+1);  typechg=double(loadframe);
        framenum=k;
        if framenum==1
            % manually place dot
            %mydot=roipoly(thisframe);
            image(loadframe);
            [mydot]=round(ginput(1));
            coords=mydot;
        else
            % pull the frame and convert
            % probably show first frame and manually ask for the led position, then
            % go from there
            
            %image(thisframe);
            %hold on;
            oldcoords=coords(framenum-skipby,:); % find the old coordinates (will be an x and y index)
            
            
            
            
            % subtract background and remove outside maze
            netimage=(typechg-meanimage).* arenamask;
            % for blue led, take max image and penalize it if its red
            % for red image, take red over other two, and add some overal luminance
            % (.2)
            highblue=netimage; highblue(:,:,3)=highblue(:,:,3)*3;
            newimage=netimage(:,:,1)-.2*mean(highblue(:,:,2:3),3);
            %imagesc(newimage>max(max(max(newimage)))*.8)
            
            
            % either zscore the whole thing, or just take a % of the strongest
            % pixel
            ratedimage=(newimage-mean(mean(newimage)))*nanstd(nanstd(newimage));
            
            % now zerp out the values outside the world of possibility
            % 50 pixel circle
            
            % place it around our rat
            closeinds= circleinds+oldcoords;
            closemask=roipoly(ratedimage,closeinds(:,1),closeinds(:,2));
            % so I think maybe i'll go from highest threshold to lower thresholds
            % until i have a single blob
            ratedimage=ratedimage.*closemask;
            % for each frame go from most stringent to least to get a good centroid
            thresh=[50:-2:1]; clear stats;
            for i=1:length(thresh)
                
                binoimage=ratedimage>thresh(i);
                %subplot(1,2,2);
                % imagesc(binoimage);
                stats=regionprops(binoimage,'area','centroid','majoraxislength','minoraxislength');
                
                % if no blobs lower threshold, if 1 blob check that its within our
                % area
                if ~exist('stats','var'), stats=[]; end
                
                if length(stats)==1
                    if stats.Area>4, break, end % got the one blob, lets use it
                elseif length(stats)>1
                    % too many blobs lets use our favorite
                    
                    clear howfar;
                    for j=1:length(stats)
                        % if its roundish and has more than 6 pixels find
                        % distance
                        if stats(j).MajorAxisLength<stats(j).MinorAxisLength*3 & stats(j).Area>4
                            deltapos=oldcoords-stats(j).Centroid;
                            howfar(j)=hypot(deltapos(1),deltapos(2));
                            
                        else % otherwise ditch this blob
                            howfar(j)=nan;
                        end
                    end
                    % now pull closest blob
                    if any(~isnan(howfar))
                        % find the closest blob (and no far away ones)
                        stats=stats(howfar==min(howfar) & howfar<20);
                        if length(stats)>1
                            
                            break
                        end
                        
                    end
                end
            end
            if ~exist('stats','var')
                fprintf('messed up on frame %d \n', framenum);
            elseif length(stats)~=1
                try
                    thresh(i)
                    [stats.Area]
                catch
                    fprintf('something really got fucked up \n');
                    keyboard
                end
                fprintf('messed up on frame %d \n', framenum);
                image(loadframe);
                drawnow;
                [mydot]=round(ginput(1));
                coords(framenum,:)=mydot;
                
            else
                coords(framenum,:)=stats.Centroid;
            end
        end
        %plot(coords(framenum,1),coords(framenum,2),'w*');
        if ~rem(framenum,(3*100)+1)
            fprintf('frame number %d done \n',framenum);
        end
    end
end
fprintf('done \n')
% then interpolate between
coords2=[(1:length(coords(:,1)))' coords];
x=coords2(sum(coords2(:,2:3),2)>0,1);
v=coords(sum(coords2(:,2:3),2)>0,:);
xq=coords2(:,1);
vq = interp1(x,v,xq,'spline');


goodcoords=[xq vq];


close
figure;
image(loadframe)
objectcoords=round(ginput(2));


coorddata.goodcoords=goodcoords;
coorddata.oldcoords=pos_data;
coorddata.objectcoords=objectcoords;

save([avipath '\' sessionname],'coorddata');

end