%Makes all .NEX files in folder readable in Matlab.
dirName=uigetdir;
fileList = getAllNEXFiles(dirName);
for i=1:length(fileList)
    disp (['Calculating ' fileList{i} ])
 
    
    [root,nex] = NEX2CMBnew(fileList{i});
    root = root.FixPos;
    save([fileList{i}(1:end-4) '_CMBobj.mat'],'root')
end