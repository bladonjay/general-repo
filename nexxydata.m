function [xydata]=nexxydata(filename)
if nargin~=1
    [filename,pathname]=uigetfile;
    filename=[pathname filename];
end

[~, ~, ~, ts, ~, coords] = nex_marker(filename, 'PLX Coord 1');
if numel(coords)<2
    [~, ~, ~, ts, ~, coords] = nex_marker(filename, 'DVT Coord 1');
end

Xcoordstr=coords(:,:,1);
Ycoordstr=coords(:,:,2);

X=nan(size(Xcoordstr,1),1);
Y=nan(size(Ycoordstr,1),1);
for m=1:size(Xcoordstr,1)
 X(m)=str2num(Xcoordstr(m,:));
 Y(m)=str2num(Ycoordstr(m,:));
end

xydata=[ts' X Y];    