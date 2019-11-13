function [outavifile] = BrightenAVI(avifile,outavifile)
%This goes through and lets you brighten your avifile, still working on the
%GUI part where you designate your brightness... 
% ***Wont write avi files with the right compression for plexon, need to
% have .avi file with JPEG2000 image compression I believe



%%%%%%%%%%%%%%
% fetch the AVI file 
if(exist('avifile','var')~=1 || isempty(avifile))
    FilterSpec = {'*.AVI', 'AVI File (*.avi)';
                  '*', 'All Files'};
    [fname, pathname] = uigetfile(FilterSpec, 'Select an input video file');
    if(fname == 0); return; end
    avifile = strcat(pathname, fname);
end

%%%%%%%%%%%%%
% create an avifile object to read from
if(ischar(avifile) && exist(avifile,'file')==2)
    vid = VideoReader(avifile);
elseif(isa(avifile,'VideoReader'))
    vid = avifile;
else
    error([mfilename ':InvalidArgument'],...
        'Forth input argument should be filename or ''VideoReader'' class.');
end
framenum=(1:vid.NumberOfFrames);


%%%%%%%%%%%%%%%%
% build the outfile to write to
dotloc = find(vid.Name=='.',1,'last');
if(exist('outavifile','var')~=1 || isempty(outavifile))
    defaultoutfile = [vid.Name(1:dotloc-1) '-brighter.avi'];
    FilterSpec = {'*.AVI', 'AVI File (*.avi)';
                  '*', 'All Files'};
    [fname, pathname] = uiputfile(FilterSpec,...
        'Select an output video file',[pathname defaultoutfile]);
    if(fname == 0); return; end
    outavifile = strcat(pathname, fname);
end

%%%%%%%%%%%%%%%%%
% imfilt filters the image, this is where you brighten, i'll have to add
% code here so the gui puts in the first frame and we can see how much we
% want to filter it
if exist('imfilt','var')
    if(~isa(imfilt,'function_handle'))
        error([mfilename ':InvalidImageFilter'],'Invalid image filter supplied.');
    end
end
%%%%%%%%%%%%%%%%%%
% run the gui to see how bright you want the video;
% pull the first frame;
testframe(:,:,:) = read(vid,framenum(1));

[finalframe,brighter]=brightengui(testframe);




%%%%%%%%%%%%%%%%%%%%
% Open the output video file. dont write anything to it
ovid = VideoWriter(outavifile,'Motion JPEG AVI');
%ovid = VideoWriter(outavifile,'MPEG-4');
open(ovid)
c1 = onCleanup(@() close(ovid));


% build waitbar
h = waitbar(0,sprintf('Please wait... (%02.0f%%)',0));
c2 = onCleanup(@(x,y) close(h));
nframes = numel(framenum);
updateintv = max(floor(nframes/100),1);


% go frame by frame
for ii = 1:nframes
    % Read the frame
    frame(:,:,:) = read(vid,framenum(ii));
    
    % Apply image filter
     frame = frame*brighter;
    writeVideo(ovid,frame);
    
    if(mod(ii,updateintv)==0)
        if(ishandle(h))
            waitbar(ii/nframes,h,sprintf('Please wait... (%02.0f%%)',ii*100/nframes));
        else break
        end
    end
end

close(ovid);
    
    
end


