[file,path]=uigetfile;


vidobj=VideoReader([path '\' file]);

%%
avgframe=double(read(vidobj,1));
ground=1:20:5000;
for i=ground;
    avgframe=avgframe+double(read(vidobj,i));
   h = waitbar(i/2000);
end
delete(h);
%%
realavg=avgframe/length(ground);


%imagesc(realavg(:,:,1))
%vidavg=uint16(realavg);
%%
figure('Position',[100,100,1800,600]);
for i=1
    
    
tryframe=double(read(vidobj,i));
subplot(1,3,1);
imagesc(tryframe(:,:,1));

subplot(1,3,2);
negframe=tryframe(:,:,1)-realavg(:,:,1);
imagesc(negframe(:,:,1),[-75 75]);
bioframe=(negframe>40) - (negframe<-10) + (negframe>60);
subplot(1,3,3);
imagesc(bioframe);
drawnow



end



%%
% algorithm;
% 1. get average image across maybe 500 discontiguous frames, maybe
% linspace across whole session
% 2. Now subtract that from each frame
% 3. for frame 1: get a boundary, then make sure rat is on maze, if not go forward five frames
% and ask again.
% 4. if yes, it will plot a center of gravity of the positive blob and one
% for the negative blob ontop of the raw image.  Ask if thats correct.  If
% not, allow to draw crossheirs.
% 5. ok save that blob, and record the coords for head and butt, and move
% to the next frame.
% 6. find the blobs in the new image, and find which ones overlap best with
% the old blob.
% 7. if there is no overlap between the current blob and the last blob,
% Look ahead one and two more frames.  If not, plot it out and use the same click gui.
% 8. keep track of the number of frames and the index of the coordinates,
% no need to attach timestamps to it, we have it in the PLX file.
%%

% first get the video

[file,path]=uigetfile;
vidobj=VideoReader([path '\' file]);
%%
% display the info
fprintf('Video name is %s \n', vidobj.name(1:end-4));
fprintf('Has %d frames and lasts %.2f minutes \n', vidobj.NumberOfFrames,vidobj.Duration/60);

% first see if image needs brightening
tryframe=read(vidobj,1);
figure(2); image(tryframe);
[brighten]=input('Does this image need brightening, y/n \n','s');
close
if brighten=='y'
    [newavi]=BrightenAVI([path '\' file], [path '\' file '-bright']);
    close(vidobj); clear vidobj;
    vidobj=VideoReader([path '\' newavi]);
end


% now plot first frame and get the mask

tryframe=read(vidobj,1);

% outline the valid portion of the maze
[xpos,ypos,tracks]=outlinemaze3(tryframe);

% now make it a polygon and get the insides
validmaze=repmat(roipoly(tryframe,xpos,ypos),[1,1,3]);
% validmaze is a mask we can use to cover the blobs in our analysis

% plot our first frame, make sure rat is on maze
figure(2); image(tryframe);
onmaze=input('Are the LED''s visible?','s');

startframe=1;
while onmaze~='y'
    figure(2); startframe=startframe+1; 
    tryframe=read(vidobj,startframe);
    image(tryframe);
    drawnow
    onmaze=input('Are the LED''s visible?','s');
end
close

%%


% now get average video
avgframe=double(read(vidobj,1));

groundframes=1:20:5000;
for i=groundframes(2:end);
    avgframe=avgframe+double(read(vidobj,i));
   h = waitbar(i/2000);
end
delete(h);
realavg=avgframe/length(groundframes);

% now start tracking
thisneg=double(tryframe)-realavg;
goodneg=thisneg.*validmaze;
% first red channel

figure('Position',[80,80,1600,400]);
subplot(1,4,1);
imagesc(goodneg(:,:,1));
thresh=prctile(linearize(goodneg(:,:,1)),99.9);
bwred=goodneg(:,:,1)>thresh;
subplot(1,4,2);
imagesc(bwred);
blobs=bwboundaries(bwred);


    subplot(1,4,3);
    newimage=zeros(size(bwred),'like',bwred);
for i=1:length(blobs)











